=head1 LICENSE

Copyright [2019] Global Gene Corp Pvt. Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 will.mclaren@globalgenecorp.com

=cut

=head1 NAME

VCF2JSON - Convert VCF input string to JSON in JSON output

=head1 SYNOPSIS

  mv VCF2JSON.pm ~/.vep/Plugins
  ./vep -i variants.vcf --plugin VCF2JSON --json

=head1 DESCRIPTION

  A VEP plugin that encodes VCF input as JSON when using
  JSON output (--json). VEP's default JSON output includes
  a string representation of the input as the value of the
  key "input" in the root of the JSON object.

  When (and only when) using VCF input (--format vep) and
  JSON output (--json), this plugin converts the value of
  the "input" key to a JSON representation of the VCF with
  the key "vcf". The value of the "input" key is deleted
  but can be retained with delete=0.

  Options, passed as e.g. --plugin VCF2JSON,opt1=val1

  - elastic=[0/1] : if set, prepends bare-bones JSON line for
    indexing record in elasticsearch (default: 1)

  - index=[index_name] : name of elasticsearch index (vep)

  - delete=[0/1] : delete the original "input" key (1)

  - info_multi=[string/list] : parse comma-separated INFO
    keys to lists or leave as strings [list]


  Notes/limitations:

  - *IMPORTANT* - fields to the right of INFO are discarded,
    this means genotype fields are *NOT* included!!!

  - Numbers are converted from strings to ints or floats
    except for the content of the CHROM field.

  - ID, ALT and FILTER fields are *always* parsed to lists
    even if they have only one member.

  - The INFO field is parsed into a sub-object/hash
    by converting the semicolon-separated key=value pairs
    to {"key": "value"}.

  - INFO key/values with multiple (comma-separated) values are
    parsed to lists. NB this may lead to inconsistencies if
    some VCF lines have only one value. Use the info_multi
    parameter (see above) to control this behaviour.

  - VCF missing values (".") are converted to JSON null

  - INFO keys with no values are given a JSON boolean
    value of true e.g. {"EXAMPLE_KEY": true}.

=cut

package VCF2JSON;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %DEFAULTS = (
  'elastic' => 1,
  'index' => 'vep',
  'info_multi' => 'list',
  'delete' => 1,
);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  my %plugin_cfg = %DEFAULTS;
  $plugin_cfg{$_} = $param_hash->{$_} for keys %$param_hash;

  $self->config->{_vcf2json_config} = \%plugin_cfg;

  return $self;
}

sub run { return {}; }

1;

package Bio::EnsEMBL::VEP::OutputFactory::JSON;

no warnings qw(redefine);

my @VCF_HEADERS = qw(CHROM POS ID REF ALT QUAL FILTER INFO);
my %MULTI = map {$_ => 1} qw(ID ALT FILTER);

sub output_hash_to_line {
  my $self = shift;
  my $hash = shift;

  $self->{json_obj} ||= JSON->new;
  my @split;

  my $cfg = $self->config->{_params}->{_vcf2json_config};

  # is input VCF?
  if(!exists($self->{_plugin_input_format})) {
    @split = split("\\t", $hash->{input});
    
    if (
      $split[0] =~ /(chr)?\w+/ &&
      $split[1] =~ /^\d+$/ &&
      $split[3] && $split[3] =~ /^[ACGTN\-\.]+$/i &&
      $split[4]
    ) {
      $self->{_plugin_input_format} = 'vcf';
    }
    else {
      $self->{_plugin_input_format} = 'other';
    }
  }

  if($self->{_plugin_input_format} eq 'vcf') {
    @split = split("\\t", $hash->{input}) unless scalar @split;

    # parse core fields (up to INFO)
    my %vcf_hash;
    for my $i(0..$#VCF_HEADERS) {
      my ($h, $v) = ($VCF_HEADERS[$i], $split[$i]);

      if($MULTI{$h}) {
        $vcf_hash{$h} = [map {_dot_to_null($_)} split(/\;|\,/, $v)];
      }
      else {
        $vcf_hash{$h} = _dot_to_null($v);
      }
    }

    # now parse INFO, key1=value1;key2=value2
    my @info_split = split(';', $vcf_hash{INFO});
    my %info;
    foreach my $i(@info_split) {
      my ($k, $v) = split('=', $i);
      
      # some keys have no value, give them a boolean value
      if(!defined($v)) {
        $info{$k} = JSON::true;
      }
      elsif($cfg->{info_multi} || '' eq 'list') {
        my @values = split(',', $v);
        if(scalar @values > 1) {
          $info{$k} = [map {_dot_to_null($_)} @values];
        }
        else {
          $info{$k} = _dot_to_null($v);
        }
      }
      else {
        $info{$k} = _dot_to_null($v);
      }
    }

    $vcf_hash{INFO} = \%info;

    # overwrite the input key in the hash sent back to VEP
    $hash->{vcf} = numberify(\%vcf_hash, {'CHROM' => 1});

    delete($hash->{input}) if $cfg->{delete};
  }

  my $ret = '';
  if($cfg->{elastic}) {
    $ret .= $self->{json_obj}->encode(
      { "index" => { "_index" => $cfg->{index}, "_type" => "_doc" } }
    )."\n";
  }
  $ret .= $self->{json_obj}->encode($hash);

  return $ret;
}

sub _dot_to_null {
  return $_[0] eq '.' ? JSON::null : $_[0];
}

1;
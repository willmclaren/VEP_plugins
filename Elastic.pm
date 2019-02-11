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

package Elastic;

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

  # enforce allele_number flag so we can map between VCF alleles and VEP output
  $self->config->{allele_number} = 1;

  return $self;
}

sub run { return {}; }

1;

package Bio::EnsEMBL::VEP::OutputFactory::JSON;
use Scalar::Util qw(looks_like_number);

no warnings qw(redefine);

my @VCF_HEADERS = qw(CHROM POS ID REF ALT QUAL FILTER INFO);
my %MULTI = map {$_ => 1} qw(ID ALT FILTER);

sub output_hash_to_line {
  my $self = shift;
  my $hash = shift;

  $self->{json_obj} ||= JSON->new;
  my $cfg = $self->config->{_params}->{_vcf2json_config};

  # is input VCF?
  if(!exists($self->{_plugin_input_format})) {
    my @split = split("\\t", $hash->{input});
    
    if (
      $split[0] =~ /(chr)?\w+/ &&
      $split[1] && $split[1] =~ /^\d+$/ &&
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
    $hash->{vcf} = $self->input_to_vcf_hash($hash->{input});
    delete($hash->{input}) if $cfg->{delete};
  }

  $hash = $self->convert_to_by_allele($hash);

  my $ret = '';
  if($cfg->{elastic}) {
    $ret .= $self->{json_obj}->encode(
      { "index" => { "_index" => $cfg->{index}, "_type" => "_doc" } }
    )."\n";
  }
  $ret .= $self->{json_obj}->encode($hash);

  return $ret;
}

sub input_to_vcf_hash {
  my ($self, $input) = @_;
  my $cfg = $self->config->{_params}->{_vcf2json_config};

  my @split = split("\\t", $input);

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
  if($vcf_hash{INFO}) {
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
  }

  return numberify(\%vcf_hash, {'CHROM' => 1});
}

sub convert_to_by_allele {
  my ($self, $hash) = @_;

  my $by_allele = {};
  my %rename = ('colocated_variants' => 'frequencies');
  my %delete = qw(transcript_consequences regulatory_feature_consequences motif_feature_consequences intergenic_consequences);

  $self->{_allele_map} = {};

  foreach my $convert_key(qw(vcf transcript_consequences regulatory_feature_consequences motif_feature_consequences intergenic_consequences colocated_variants)) {
    my $method_name = 'convert_'.$convert_key;
    my $converted = $self->$method_name($hash->{$convert_key} || []) if $self->can($method_name);

    $by_allele->{$_}->{$rename{$convert_key} || $convert_key} = $converted->{$_} for keys %$converted;

    delete $hash->{$convert_key} if $delete{$convert_key};
  }

  # convert string allele keys to index
  foreach my $a(keys %$by_allele) {
    if(!looks_like_number($a)) {
      if(my $i = $self->{_allele_map}->{$a}) {
        my $d = delete $by_allele->{$a};

        $by_allele->{$i}->{$_} = $d->{$_} for keys %$d;
      }
    }
  }

  my @alts = @{$self->{_alts}};
  my %reverse_map = map {$self->{_allele_map}->{$_} => $_} keys %{$self->{_allele_map}};
  for my $i(0..$#alts) {
    my $data = {
      'alt' => $alts[$i],
      'vep_allele' => $reverse_map{$i + 1} || '?',
      'allele_num' => $i + 1
    };

    $data->{$_} = $by_allele->{$i + 1}->{$_} for keys %{$by_allele->{$i + 1} || {}};
    push @{$hash->{by_allele}}, numberify($data);
  }

  return $hash;
}

sub convert_vcf {
  my ($self, $data) = @_;
  return {} if ref($data) eq 'ARRAY' && scalar(@$data) == 0;

  my @alts = @{$data->{ALT}};
  my $by_allele;

  for my $key(qw(AN AF AC AC_Het AC_Hom AC_Hemi HWE ExcHet)) {
    if(my $val = delete $data->{INFO}->{$key}) {
      if(ref($val) ne 'ARRAY') {
        $val = [$val];
      }
      $by_allele->{$_ + 1}->{$key} = $val->[$_] for 0..$#alts;
    }
  }

  $self->{_alts} = $data->{ALT};

  return $by_allele;
}

sub convert_colocated_variants {
  my ($self, $data) = @_;

  my $by_allele;

  foreach my $v(@$data) {
    if(my $freqs = $v->{frequencies}) {
      $by_allele = $freqs;
      delete $v->{frequencies};
    }
  }

  return $by_allele;
}

sub convert_transcript_consequences {
  my $self = shift;
  return $self->_generic_convert_consequences(@_);
}

sub convert_regulatory_feature_consequences {
  my $self = shift;
  return $self->_generic_convert_consequences(@_);
}

sub convert_motif_feature_consequences {
  my $self = shift;
  return $self->_generic_convert_consequences(@_);
}

sub convert_intergenic_consequences {
  my $self = shift;
  return $self->_generic_convert_consequences(@_);
}

sub _generic_convert_consequences {
  my ($self, $data) = @_;

  my $by_allele;

  foreach my $tc(@$data) {
    my $num = delete $tc->{allele_num};
    $self->{_allele_map}->{delete $tc->{variant_allele}} = $num;

    push @{$by_allele->{$num}}, $tc;
  }

  return $by_allele;
}

sub _dot_to_null {
  return $_[0] eq '.' ? JSON::null : $_[0];
}

1;
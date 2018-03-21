#!/usr/bin/env bats

@test "SRA Toolkit v2.8.2" {
  v="$(fastq-dump --version)"
  [[ "$v" =~ "2.8.2" ]]
}


@test "AWS CLI v1.11.146" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "1.11.146" ]]
}


@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
}

@test "MetaSPAdes v3.11.1" {
  v="$(metaspades.py -h 2>&1)"

  [[ "$v" =~ "v3.11.1" ]]
}

@test "Prokka v1.12" {
  v="$(prokka --version 2>&1)"

  [[ "$v" =~ "1.12" ]]
}

@test "tbl2asn" {

  [[ "$( tbl2asn 2>&1 || true )" =~ "[tbl2asn]" ]]

}

@test "Make sure the run script is in the PATH" {
  h="$(run.py -h)"

  [[ "$h" =~ "Assemble a set of reads using metaSPAdes" ]]
}

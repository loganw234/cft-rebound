# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Prove what a packaged kernel IP instantiates, without linking it.
# Instantiate the packaged IP in a throwaway project exactly the way the
# Vitis link project will, generate its synthesis wrapper, and print the
# parameter overrides that wrapper hands to cft_krnl. A generic that did
# not survive packaging shows up here as its RTL default, twenty minutes
# after packaging rather than two hours after a link.
#
#   vivado -mode batch -source verify_variant_xo.tcl -tclargs <part> <pkg_dir> <work_dir>

set part     "xcu50-fsvh2104-2-e"
set pkg_dir  ""
set work_dir ""
if {$argc >= 1} { set part     [lindex $argv 0] }
if {$argc >= 2} { set pkg_dir  [file normalize [lindex $argv 1]] }
if {$argc >= 3} { set work_dir [file normalize [lindex $argv 2]] }
if {$pkg_dir eq "" || $work_dir eq ""} {
  error "usage: -tclargs <part> <pkg_dir> <work_dir>"
}

create_project -force xo_verify $work_dir -part $part
set_property ip_repo_paths $pkg_dir [current_project]
update_ip_catalog
create_ip -vlnv improperaperture.com:kernel:cft_krnl:1.0 -module_name cft_krnl_v
generate_target {instantiation_template synthesis} [get_ips cft_krnl_v]

# The wrapper the IP integrator would synthesise.
set found 0
foreach f [get_files -all -of_objects [get_ips cft_krnl_v]] {
  if {[string match "*/synth/cft_krnl_v.v" $f] || [string match "*/synth/cft_krnl_v.sv" $f]} {
    set found 1
    puts "WRAPPER: $f"
    set fh [open $f r]
    set txt [read $fh]
    close $fh
    foreach line [split $txt "\n"] {
      if {[regexp {\.(EN_FP64|EN_FP128|EN_FP256|BEAT_BITS|MUL_PASSES|FUSE_[A-Z]+|AR_DEPTH|AW_DEPTH|FIFO_LOG2|BURST_LOG2)\s*\(} $line]} {
        puts "WRAPPER_PARAM: [string trim $line]"
      }
    }
  }
}
if {!$found} { error "no synthesis wrapper generated for cft_krnl_v" }
close_project

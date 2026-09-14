# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Package cft-fp256's cft_krnl into a Vitis kernel object (.xo) with
# GENERIC OVERRIDES - the one thing cft-fp256's own hw/package_kernel.tcl
# cannot do. That script removes every user parameter from the packaged
# IP, so a bitstream can only ever carry the RTL defaults (the full
# fp32/fp64/fp128/fp256 tile). This one sets the requested values on the
# HDL parameters BEFORE the user parameters are removed, so the value the
# IP integrator instantiates is the override, and prints every HDL
# parameter it packaged so the log says what the image carries.
#
#   CFT_GENERICS="EN_FP256=0" vivado -mode batch -source package_variant.tcl \
#       -tclargs <part> <build_dir> <rtl_dir> <kernel_xml>
#
# Generics come from the environment, not -tclargs: on Windows vivado.bat
# is a cmd.exe wrapper and cmd.exe splits an argument at `=`, which
# cft-fp256's hw/impl_krnl_ooc.tcl learned the hard way (a parameter with
# no value, bound to 0, 2026-09-06). Kernel name and VLNV stay cft_krnl:
# libcft's XRT backend opens compute units by the literal name
# `cft_krnl:{cft_krnl_N}` (host/src/backend_xrt.cpp), so a renamed
# variant would open zero tiles. What the image lacks is published by
# the tile itself, in CAPS[3:0], and refused by the tile with STATUS[3].
#
# Everything else follows package_kernel.tcl step for step (the
# Vitis_Accel_Examples RTL kernel sequence), so the only difference
# between this .xo and cft-fp256's is the parameter values.

set part       "xcu50-fsvh2104-2-e"
set build_dir  "build"
set rtl_dir    ""
set kernel_xml ""
if {$argc >= 1} { set part       [lindex $argv 0] }
if {$argc >= 2} { set build_dir  [lindex $argv 1] }
if {$argc >= 3} { set rtl_dir    [lindex $argv 2] }
if {$argc >= 4} { set kernel_xml [lindex $argv 3] }
if {$rtl_dir eq "" || $kernel_xml eq ""} {
  error "usage: -tclargs <part> <build_dir> <rtl_dir> <kernel_xml>  (generics in CFT_GENERICS)"
}
set generics {}
if {[info exists ::env(CFT_GENERICS)]} { set generics $::env(CFT_GENERICS) }
if {[llength $generics] == 0} {
  # Refuse rather than quietly package the full tile under a variant's
  # name: an empty override list is exactly the mistake this file exists
  # to make impossible.
  error "CFT_GENERICS is empty; to package the full tile use cft-fp256's hw/package_kernel.tcl"
}

set rtl_dir [file normalize $rtl_dir]
set pkg_dir [file normalize "$build_dir/packaged_kernel"]
set tmp_dir [file normalize "$build_dir/tmp_kernel_pack"]
set xo_path [file normalize "$build_dir/cft_krnl.xo"]
puts "VARIANT: part=$part rtl=$rtl_dir generics={$generics}"

file mkdir $build_dir

create_project -force kernel_pack $tmp_dir -part $part
add_files -norecurse [glob $rtl_dir/*.sv]
set_property top cft_krnl [current_fileset]
update_compile_order -fileset sources_1

ipx::package_project -root_dir $pkg_dir -vendor improperaperture.com \
    -library kernel -taxonomy /KernelIP -import_files -set_current false
ipx::unload_core $pkg_dir/component.xml
ipx::edit_ip_in_project -upgrade true -name tmp_edit_project \
    -directory $pkg_dir $pkg_dir/component.xml

set core [ipx::current_core]
set_property core_revision 2 $core

# ---- the overrides ----------------------------------------------------
# Set on the HDL parameter (what the instantiation carries) and on the
# user parameter while it still exists (so the two cannot disagree in
# the moment between here and the removal below). A name the RTL does
# not declare is an error, not a warning: a typo would otherwise package
# the full tile and say nothing.
foreach g $generics {
  if {![regexp {^([A-Za-z_][A-Za-z0-9_]*)=(.+)$} $g -> name value]} {
    error "generic '$g' is not NAME=VALUE"
  }
  set hp [ipx::get_hdl_parameters $name -of_objects $core]
  if {$hp eq ""} { error "cft_krnl declares no parameter named $name" }
  set old [get_property value $hp]
  # A `parameter bit` is stored as a QUOTED bit string - the packager
  # writes EN_FP64's default as "1" (format bitString, length 1), and
  # the 2026.1 packager accepted a bare 0 too but then carried two
  # spellings for one type. Keep the packager's own spelling: quoted,
  # same length, digits substituted. Anything else (int, long) is taken
  # verbatim.
  if {[regexp {^"([01]+)"$} $old -> bits]} {
    if {![regexp {^[01]+$} $value] || [string length $value] != [string length $bits]} {
      error "$name is a ${[string length $bits]}-bit parameter; '$value' is not a value for it"
    }
    set value "\"$value\""
  } elseif {[regexp {^(\d+)'b[01]+$} $old -> w]} {
    set value "${w}'b${value}"
  }
  set up [ipx::get_user_parameters $name -of_objects $core]
  if {$up ne ""} { set_property value $value $up }
  set_property value $value $hp
  puts "GENERIC: $name  $old -> [get_property value $hp]"
}

set_property sdx_kernel true $core
set_property sdx_kernel_type rtl $core
set_property supported_families { } $core
set_property auto_family_support_level level_2 $core
foreach up [ipx::get_user_parameters -of_objects $core] {
  ipx::remove_user_parameter [get_property NAME $up] $core
}
# What the .xo carries, by name, after the removal - this is the line a
# manifest and a reader should trust, not the intent above.
foreach hp [ipx::get_hdl_parameters -of_objects $core] {
  puts "HDLPARAM: [get_property NAME $hp] = [get_property value $hp]"
}
ipx::create_xgui_files $core
ipx::associate_bus_interfaces -busif s_axi_control -clock ap_clk $core
foreach busif {m_axi_a m_axi_b m_axi_c m_axi_d} {
  ipx::associate_bus_interfaces -busif $busif -clock ap_clk $core
}
ipx::update_checksums $core
ipx::save_core $core
close_project -delete

package_xo -xo_path $xo_path -kernel_name cft_krnl \
    -ip_directory $pkg_dir -kernel_xml $kernel_xml -force

puts "INFO: wrote $xo_path"

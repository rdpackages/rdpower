version 16.0
clear all
set more off

args repo_root output
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
global RDP_REPO_ROOT "`repo_root'"
if `"`output'"' == "" {
    local output "`repo_root'/docs/audit/baselines/stata-current.json"
}
local output : subinstr local output "\" "/", all

capture mkdir "`repo_root'/docs"
capture mkdir "`repo_root'/docs/audit"
capture mkdir "`repo_root'/docs/audit/baselines"

adopath ++ "`repo_root'/stata"
capture confirm file "`repo_root'/../rdrobust/stata/rdrobust.ado"
if !_rc {
    adopath ++ "`repo_root'/../rdrobust/stata"
}

program define preload_mata
    capture confirm file "$RDP_REPO_ROOT/../rdrobust/stata/rdrobust_functions.do"
    if !_rc {
        quietly do "$RDP_REPO_ROOT/../rdrobust/stata/rdrobust_functions.do"
    }
end

program define write_rdpower_case
    args handle casename comma include_conv
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "rdpower": {"' _n
    file write `handle' `"        "power_rbc": "' %24.16e (r(power_rbc)) `","' _n
    file write `handle' `"        "se_rbc": "' %24.16e (r(se_rbc)) `","' _n
    if "`include_conv'" == "conv" {
        file write `handle' `"        "power_conv": "' %24.16e (r(power_conv)) `","' _n
        file write `handle' `"        "se_conv": "' %24.16e (r(se_conv)) `","' _n
    }
    file write `handle' `"        "sampsi_l": "' %24.16e (r(sampsi_l)) `","' _n
    file write `handle' `"        "sampsi_r": "' %24.16e (r(sampsi_r)) `","' _n
    file write `handle' `"        "samph_l": "' %24.16e (r(samph_l)) `","' _n
    file write `handle' `"        "samph_r": "' %24.16e (r(samph_r)) `","' _n
    file write `handle' `"        "N_l": "' %24.16e (r(N_l)) `","' _n
    file write `handle' `"        "N_r": "' %24.16e (r(N_r)) `","' _n
    file write `handle' `"        "Nh_l": "' %24.16e (r(N_h_l)) `","' _n
    file write `handle' `"        "Nh_r": "' %24.16e (r(N_h_r)) `","' _n
    file write `handle' `"        "tau": "' %24.16e (r(tau)) `","' _n
    file write `handle' `"        "bias_l": "' %24.16e (r(bias_l)) `","' _n
    file write `handle' `"        "bias_r": "' %24.16e (r(bias_r)) `","' _n
    file write `handle' `"        "Vl_rb": "' %24.16e (r(Vl_rb)) `","' _n
    file write `handle' `"        "Vr_rb": "' %24.16e (r(Vr_rb)) `","' _n
    file write `handle' `"        "alpha": "' %24.16e (r(alpha)) _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

program define write_rdsampsi_case
    args handle casename comma include_conv
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "rdsampsi": {"' _n
    file write `handle' `"        "sampsi.h.tot": "' %24.16e (r(sampsi_h_tot)) `","' _n
    file write `handle' `"        "sampsi.h.l": "' %24.16e (r(sampsi_h_l)) `","' _n
    file write `handle' `"        "sampsi.h.r": "' %24.16e (r(sampsi_h_r)) `","' _n
    file write `handle' `"        "sampsi.tot": "' %24.16e (r(sampsi_tot)) `","' _n
    if "`include_conv'" == "conv" {
        file write `handle' `"        "sampsi.h.tot.cl": "' %24.16e (r(sampsi_h_tot_cl)) `","' _n
        file write `handle' `"        "sampsi.h.l.cl": "' %24.16e (r(sampsi_h_l_cl)) `","' _n
        file write `handle' `"        "sampsi.h.r.cl": "' %24.16e (r(sampsi_h_r_cl)) `","' _n
        file write `handle' `"        "sampsi.tot.cl": "' %24.16e (r(sampsi_tot_cl)) `","' _n
    }
    file write `handle' `"        "N.l": "' %24.16e (r(N_l)) `","' _n
    file write `handle' `"        "N.r": "' %24.16e (r(N_r)) `","' _n
    file write `handle' `"        "Nh.l": "' %24.16e (r(N_h_l)) `","' _n
    file write `handle' `"        "Nh.r": "' %24.16e (r(N_h_r)) `","' _n
    file write `handle' `"        "bias.l": "' %24.16e (r(bias_l)) `","' _n
    file write `handle' `"        "bias.r": "' %24.16e (r(bias_r)) `","' _n
    file write `handle' `"        "var.l": "' %24.16e (r(var_l)) `","' _n
    file write `handle' `"        "var.r": "' %24.16e (r(var_r)) `","' _n
    file write `handle' `"        "samph.l": "' %24.16e (r(samph_l)) `","' _n
    file write `handle' `"        "samph.r": "' %24.16e (r(samph_r)) `","' _n
    file write `handle' `"        "tau": "' %24.16e (r(tau)) `","' _n
    file write `handle' `"        "beta": "' %24.16e (r(beta)) `","' _n
    file write `handle' `"        "alpha": "' %24.16e (r(alpha)) `","' _n
    file write `handle' `"        "init.cond": "' %24.16e (r(init_cond)) `","' _n
    file write `handle' `"        "no.iter": "' %24.16e (r(no_iter)) _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

program define write_rdmde_case
    args handle casename comma
    file write `handle' `"`comma'    "`casename'": {"' _n
    file write `handle' `"      "rdmde": {"' _n
    file write `handle' `"        "mde": "' %24.16e (r(mde)) `","' _n
    file write `handle' `"        "se_rbc": "' %24.16e (r(se_rbc)) `","' _n
    file write `handle' `"        "sampsi_l": "' %24.16e (r(sampsi_l)) `","' _n
    file write `handle' `"        "sampsi_r": "' %24.16e (r(sampsi_r)) `","' _n
    file write `handle' `"        "samph_l": "' %24.16e (r(samph_l)) `","' _n
    file write `handle' `"        "samph_r": "' %24.16e (r(samph_r)) `","' _n
    file write `handle' `"        "N_l": "' %24.16e (r(N_l)) `","' _n
    file write `handle' `"        "N_r": "' %24.16e (r(N_r)) `","' _n
    file write `handle' `"        "Nh_l": "' %24.16e (r(N_h_l)) `","' _n
    file write `handle' `"        "Nh_r": "' %24.16e (r(N_h_r)) `","' _n
    file write `handle' `"        "bias_l": "' %24.16e (r(bias_l)) `","' _n
    file write `handle' `"        "bias_r": "' %24.16e (r(bias_r)) `","' _n
    file write `handle' `"        "Vl_rb": "' %24.16e (r(Vl_rb)) `","' _n
    file write `handle' `"        "Vr_rb": "' %24.16e (r(Vr_rb)) `","' _n
    file write `handle' `"        "alpha": "' %24.16e (r(alpha)) `","' _n
    file write `handle' `"        "beta": "' %24.16e (r(beta)) _n
    file write `handle' `"      }"' _n
    file write `handle' `"    }"' _n
end

file open jout using "`output'", write replace
file write jout `"{"' _n
file write jout `"  "schema_version": 1,"' _n
file write jout `"  "package": "rdpower","' _n
file write jout `"  "language": "stata","' _n
file write jout `"  "source": "working-tree","' _n
file write jout `"  "timestamp_utc": null,"' _n
file write jout `"  "environment": {"stata_version": "`c(stata_version)'", "platform": "`c(os)'"}, "' _n
file write jout `"  "cases": {"' _n

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdpow demvoteshfor2 demmv, tau(5)
write_rdpower_case jout senate_default ""

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdsampsi demvoteshfor2 demmv, tau(5)
write_rdsampsi_case jout senate_default_sampsi ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdmde demvoteshfor2 demmv
write_rdmde_case jout senate_default_mde ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdpow demvoteshfor2 demmv, tau(5) covs(population dopen dmidterm)
write_rdpower_case jout senate_covs ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdsampsi demvoteshfor2 demmv, tau(5) covs(population dopen dmidterm)
write_rdsampsi_case jout senate_covs_sampsi ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdmde demvoteshfor2 demmv, covs(population dopen dmidterm)
write_rdmde_case jout senate_covs_mde ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdpow demvoteshfor2 demmv, tau(5) h(16 18) b(18 20) all
write_rdpower_case jout senate_fixed_bandwidth "," conv

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdsampsi demvoteshfor2 demmv, tau(5) beta(.9) samph(18 19) nratio(.5) all
write_rdsampsi_case jout senate_fixed_bandwidth_sampsi "," conv

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdmde demvoteshfor2 demmv, beta(.75) samph(12 13)
write_rdmde_case jout senate_fixed_bandwidth_mde ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdpow demvoteshfor2 demmv, tau(5) vce(cluster state)
write_rdpower_case jout senate_cluster ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdsampsi demvoteshfor2 demmv, tau(5) vce(cluster state)
write_rdsampsi_case jout senate_cluster_sampsi ","

use "`repo_root'/stata/rdpower_senate.dta", clear
preload_mata
rdmde demvoteshfor2 demmv, vce(cluster state)
write_rdmde_case jout senate_cluster_mde ","

file write jout `"  }"' _n
file write jout `"}"' _n
file close jout

display as text "Wrote `output'"

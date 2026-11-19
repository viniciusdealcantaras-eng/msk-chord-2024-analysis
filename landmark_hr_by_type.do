* Landmark com entrada tardia e HR do estadio IV por tipo tumoral (MSK-CHORD)
* Entrada: out/analytic_dataset.csv gerado por make_analytic_dataset.py
* Nao salva nem sobrescreve nenhum banco.

version 16
clear all
set more off

import delimited "out/analytic_dataset.csv", clear varnames(1)
keep if included == 1

gen double T  = os_months_from_diagnosis
gen double L  = dx_to_seq_months
gen byte   E  = os_status
gen byte stage4 = stage_highest_recorded == "IV"
encode cancer_type, gen(ct)

* ---- Modelo principal com truncamento a esquerda (Tabela 2) ----
stset T, failure(E == 1) enter(time L)
stcox stage4, strata(ct)

* ---- HR do estadio IV por tipo tumoral (Tabela S4) ----
levelsof ct, local(tipos)
foreach k of local tipos {
    display as text _n "Tipo: `: label (ct) `k''"
    stcox stage4 if ct == `k'
}

* ---- Landmark 6 e 12 meses, entrada tardia para quem foi sequenciado depois (Tabela S3) ----
foreach t0 in 6 12 {
    preserve
    keep if T > `t0'
    gen double en = max(L, `t0')
    stset T, failure(E == 1) enter(time en)
    display as text _n "Landmark `t0' meses"
    count
    stci, median
    stcox stage4, strata(ct)
    restore
}

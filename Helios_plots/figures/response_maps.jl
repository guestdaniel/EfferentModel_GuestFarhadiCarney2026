using CairoMakie
using AuditorySignalUtils
using Helios_plots
using Helios

m = AuditorySubcortexGFCv1(
    ; 
    cf=[2000.0], 
    stage="hsr", 
    fractional=false,
    moc_weight_wdr=1.0,
    moc_weight_ic=1.0,
)
a = AnalysisRMModel(
    m; 
    freqs=LogRange(2000.0 * 2^-2, 2000.0 * 2^2, 40),
    levels=[10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0],
    dur=0.1,
)
load_or_run!(a)
viz(a.rm)
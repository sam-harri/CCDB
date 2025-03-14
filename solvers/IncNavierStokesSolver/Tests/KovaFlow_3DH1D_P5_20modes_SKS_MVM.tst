<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.11934e-06</value>
            <value variable="v" tolerance="1e-11">1.39541e-06</value>
            <value variable="w" tolerance="1e-11">1.31885e-06</value>
	    <value variable="p" tolerance="1e-10">2.9429e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">2.87017e-06</value>
            <value variable="v" tolerance="1e-11">2.25226e-06</value>
            <value variable="w" tolerance="1e-11">1.95419e-06</value>
	    <value variable="p" tolerance="1e-10">6.85396e-05</value>
        </metric>
    </metrics>
</test>

<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_Deal.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_Deal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.20841e-06</value>
            <value variable="v" tolerance="1e-11">1.33289e-06</value>
            <value variable="w" tolerance="1e-11">1.27772e-06</value>
	    <value variable="p" tolerance="1e-10">2.87963e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">3.14611e-06</value>
            <value variable="v" tolerance="1e-11">2.25323e-06</value>
            <value variable="w" tolerance="1e-11">1.80613e-06</value>
	    <value variable="p" tolerance="1e-11">6.84369e-05</value>
        </metric>
    </metrics>
</test>

<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">2.18753e-06</value>
            <value variable="v" tolerance="1e-11">1.34437e-06</value>
            <value variable="w" tolerance="1e-11">1.27768e-06</value>
	    <value variable="p" tolerance="1e-10">2.89049e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">3.14458e-06</value>
            <value variable="v" tolerance="1e-11">2.25387e-06</value>
            <value variable="w" tolerance="1e-11">1.8079e-06</value>
	    <value variable="p" tolerance="1e-10">6.87022e-05</value>
        </metric>
    </metrics>
</test>

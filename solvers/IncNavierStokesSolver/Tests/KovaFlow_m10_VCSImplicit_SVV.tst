<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay flow test case for using SVV with VCSImplicit</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_VCSImplicit.xml -I SpectralVanishingViscosity=DGKernel</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_VCSImplicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-09">7.50800e-07</value>
            <value variable="v" tolerance="1e-09">4.30729e-07</value>
            <value variable="p" tolerance="1e-09">2.63507e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.60404e-06</value>
            <value variable="v" tolerance="1e-09">7.98766e-07</value>
            <value variable="p" tolerance="1e-08">2.65873e-06</value>
        </metric>
    </metrics>
</test>

<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of a self-propelled NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>SelfPropel_MRF_naca0012.xml</parameters>
    <files>
        <file description="Session File">SelfPropel_MRF_naca0012.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">0.357952</value>
            <value variable="v" tolerance="5e-4">0.405874</value>
            <value variable="p" tolerance="5e-2">0.676777</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-3">1.58646</value>
            <value variable="v" tolerance="5e-4">1.46267</value>
            <value variable="p" tolerance="5e-2">5.60387</value>
        </metric>
    </metrics>
</test>

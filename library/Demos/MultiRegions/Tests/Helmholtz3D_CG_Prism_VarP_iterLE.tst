<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for Prism</description>
    <executable>Helmholtz3D</executable>
    <parameters>-I GlobalSysSoln=IterativeStaticCond -I Preconditioner=LowEnergyBlock Helmholtz3D_Prism.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Prism.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">0.000385246</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">0.00539787</value>
        </metric>
    </metrics>
</test>



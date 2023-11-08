<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using (Semi-Implicit) GJP Stabilisation and dealiasing with the Implicit VCS</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_GJP.xml -I SolverType=VCSImplicit</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-4">0.00111739</value>
            <value variable="v" tolerance="1e-5">0.000229595</value>
            <value variable="w" tolerance="1e-5">0.000213419</value>
            <value variable="p" tolerance="1e-4">0.00129578</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-4">0.00900795</value>
            <value variable="v" tolerance="1e-4">0.0014178</value>
            <value variable="w" tolerance="1e-4">0.00105719</value>
            <value variable="p" tolerance="1e-4">0.00806811</value>
        </metric>
    </metrics>
</test>

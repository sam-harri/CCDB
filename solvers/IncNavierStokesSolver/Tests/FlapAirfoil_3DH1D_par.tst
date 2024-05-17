<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with flapping NACA0012 airfoil</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>NACA0012_Re400.xml FlapAirfoil_3DH1D_par.xml --npz 2</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">NACA0012_Re400.xml</file>
        <file description="Session File">FlapAirfoil_3DH1D_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">31.6227</value>
            <value variable="v" tolerance="1e-6">0.137227</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">9.2039</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">3.0018</value>
            <value variable="v" tolerance="1e-6">2.95194</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">11.8778</value>
        </metric>
    </metrics>
</test>

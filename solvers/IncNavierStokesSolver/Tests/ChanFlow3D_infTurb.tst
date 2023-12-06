<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fully 3D Channel Flow. Synthetic turbulence generation is introduced in the flow field using forcing tag. </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow3D_infTurb.xml</parameters>
    <files>
        <file description="Session File">ChanFlow3D_infTurb.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">3.64494</value>
            <value variable="v" tolerance="1e-8">0.00347134</value>
            <value variable="w" tolerance="1e-6">0.001361</value>
            <value variable="p" tolerance="1e-7">0.0183825</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">1.19138</value>
            <value variable="v" tolerance="1e-7">0.0103911</value>
            <value variable="w" tolerance="1e-8">0.00984826</value>
            <value variable="p" tolerance="1e-7">0.0168057</value>
        </metric>
    </metrics>
</test>



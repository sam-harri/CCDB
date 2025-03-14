<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D stretch output </description>
    <executable python="true">chan3DH1D_stretch.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
        <file description="Session File">chan3DH1D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.816497</value>
            <value variable="y" tolerance="1e-6">0.816497</value>
            <value variable="z" tolerance="1e-5">1.55624</value>
            <value variable="u" tolerance="1e-6">0.258199</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.63299</value>
        </metric>
    </metrics>
</test>


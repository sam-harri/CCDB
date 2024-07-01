<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process Moving Reference Frame field </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m MRF:vectors=u,v MRF.xml MRF.fld MRF.plt</parameters>
    <files>
        <file description="Session File">MRF.xml</file>
        <file description="Session File">MRF.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x"   tolerance="1e-4">61.2372</value>
            <value variable="y"   tolerance="1e-4">35.3553</value>
            <value variable="u"   tolerance="1e-4">10.6066</value>
            <value variable="v" tolerance="1e-4">6.12372</value>
            <value variable="p" tolerance="1e-4">0</value>
            <value variable="xCoord" tolerance="1e-4">55.9017</value>
            <value variable="yCoord" tolerance="1e-4">43.3013</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"   tolerance="1e-4">10</value>
            <value variable="y"   tolerance="1e-4">5</value>
            <value variable="u"   tolerance="1e-4">0.866033</value>
            <value variable="v" tolerance="1e-4">0.500005</value>
            <value variable="p" tolerance="1e-4">0</value>
            <value variable="xCoord" tolerance="1e-4">11.1603</value>
            <value variable="yCoord" tolerance="1e-4">9.33013</value>
        </metric>
    </metrics>
</test>


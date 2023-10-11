<?xml version="1.0" encoding="utf-8"?>
<test runs="20">
    <description>MMF SWE solver, DG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>Perf_TestMMFSWEPlane.xml</parameters>
    <files>
        <file description="Session File">Perf_TestMMFSWEPlane.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-8">0.00109671</value>
            <value variable="u" tolerance="1e-8">0.00123931</value>
            <value variable="v" tolerance="1e-8">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-5">0.032063</value>
            <value variable="u" tolerance="1e-5">0.041687</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="1e-4" hostname="42.debian-bullseye-performance-build-and-test">0.00210637</value>
        </metric>
    </metrics>
</test>



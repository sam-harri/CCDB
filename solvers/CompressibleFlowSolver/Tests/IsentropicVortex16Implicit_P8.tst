<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16Implicit_P8.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16Implicit_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-07">9.08025e-07</value>
            <value variable="rhou" tolerance="1e-07">2.14437e-06</value>
            <value variable="rhov" tolerance="1e-07">2.09341e-06</value>
            <value variable="E" tolerance="1e-07">7.16325e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-07">5.72323e-06</value>
            <value variable="rhou" tolerance="1e-07">6.91781e-06</value>
            <value variable="rhov" tolerance="1e-07">8.37104e-06</value>
            <value variable="E" tolerance="1e-06">2.41107e-05</value>
        </metric>
    </metrics>
</test>

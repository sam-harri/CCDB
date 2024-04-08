<?xml version="1.0" encoding="utf-8"?>
<test runs="10">
    <description>Euler Isentropic Vortex P=8, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Perf_IsentropicVortex16Implicit_P8.xml</parameters>
    <files>
        <file description="Session File">Perf_IsentropicVortex16Implicit_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">2.81154e-06</value>
            <value variable="rhou" tolerance="1e-12">4.45389e-06</value>
            <value variable="rhov" tolerance="1e-12">6.25694e-06</value>
            <value variable="E" tolerance="1e-12">1.4556e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.04732e-05</value>
            <value variable="rhou" tolerance="1e-12">1.0708e-05</value>
            <value variable="rhov" tolerance="1e-12">9.07312e-06</value>
            <value variable="E" tolerance="1e-12">3.03058e-05</value>
        </metric>
        <metric type="ExecutionTime" id="3">
            <value tolerance="1e0" hostname="42.debian-bullseye-performance-build-and-test">36.0</value>
        </metric>
    </metrics>
</test>

<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady FRDG inviscidBurgers GLL_LAGRANGE_SEM, P=10</description>
    <executable>ADRSolver</executable>
    <parameters>InviscidBurgers1D_FRDG_GLL_LAGRANGE_SEM.xml</parameters>
    <files>
        <file description="Session File">InviscidBurgers1D_FRDG_GLL_LAGRANGE_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-5">2.28217</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-1">2.49987</value>
        </metric>
    </metrics>
</test>

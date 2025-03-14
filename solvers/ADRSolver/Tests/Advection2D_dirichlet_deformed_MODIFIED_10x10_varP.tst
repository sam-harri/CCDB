<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady WeakDG advection MODIFIED, P=3, Dirichlet bcs, deformed elements, variable orders</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_dirichlet_deformed_MODIFIED_10x10_varP.xml</parameters>
    <files>
        <file description="Session File">Advection2D_dirichlet_deformed_MODIFIED_10x10_varP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.00010939 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.00111027 </value>
        </metric>
    </metrics>
</test>

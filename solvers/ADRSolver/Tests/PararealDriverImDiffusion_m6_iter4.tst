<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 1, P=6 </description>
    <executable>ADRSolver</executable>
    <parameters>--npt 8 PararealDriverImDiffusion_m6_iter4.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> PararealDriverImDiffusion_m6_iter4.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-09">0.000127434</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-09">0.000304225</value>
        </metric>
    </metrics>
</test>

<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Advection FRDG Diffusion LFRDG MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_FRDG_LFRDG_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_FRDG_LFRDG_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-13">5.25674e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.06493e-07</value>
        </metric>
    </metrics>
</test>

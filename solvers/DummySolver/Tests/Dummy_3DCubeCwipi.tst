<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Run two CWIPI-coupled instances of the DummySolver, exchanging fields</description>
    <segment type="parallel">
        <executable> DummySolver </executable>
        <parameters> --verbose --cwipi 'Dummy0' Dummy_3DCubeCwipi_0.xml cube.left.xml </parameters>
        <processes> 1 </processes>
    </segment>
    <segment type="parallel">
        <executable> DummySolver </executable>
        <parameters> --verbose --cwipi 'Dummy1' Dummy_3DCubeCwipi_1.xml cube.right.xml </parameters>
        <processes> 1 </processes>
    </segment>
    <files>
        <file description="Session File 0">Dummy_3DCubeCwipi_0.xml</file>
        <file description="Session File 1">Dummy_3DCubeCwipi_1.xml</file>
        <file description="Mesh File">cube.left.xml</file>
        <file description="Mesh File">cube.right.xml</file>
        <file description="PETSc config file">.petscrc</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>.*L 2 error\s*\(variable\s*F_0_\w0\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="F_0_u0S" tolerance="5E-2">10.5</field>
                </match>
                <match>
                    <field variable="F_0_v0S" tolerance="5E-2">15.5</field>
                </match>
                <match>
                    <field variable="F_0_u0R" tolerance="5E-2">30.5</field>
                </match>
                <match>
                    <field variable="F_0_v0R" tolerance="5E-2">35.5</field>
                </match>

            </matches>
        </metric>
        <metric type="regex" id="1">
            <regex>.*L 2 error\s*\(variable\s*F_0_\w1\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="F_0_u1S" tolerance="5E-2">30.50</field>
                </match>
                <match>
                    <field variable="F_0_v1S" tolerance="5E-2">35.50</field>
                </match>
                <match>
                    <field variable="F_0_u1R" tolerance="5E-3">22.1497</field>
                </match>
                <match>
                    <field variable="F_0_v1R" tolerance="5E-3">26.4196</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>.*L 2 error\s*\(variable\s*\w0\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u0S" tolerance="5E-7">0.000122474</field>
                </match>
                <match>
                    <field variable="v0S" tolerance="5E-7">0.000122474</field>
                </match>
                <match>
                    <field variable="u0R" tolerance="5E-6">1.03161e-05</field>
                </match>
                <match>
                    <field variable="v0R" tolerance="5E-6">1.03161e-05</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="3">
            <regex>.*L 2 error\s*\(variable\s*\w1\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u1S" tolerance="5E-7">1.2247e-05</field>
                </match>
                <match>
                    <field variable="v1S" tolerance="5E-7">1.2247e-05</field>
                </match>
                <match>
                    <field variable="u1R" tolerance="5E-8">5.15227e-05</field>
                </match>
                <match>
                    <field variable="v1R" tolerance="5E-8">5.15227e-05</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>


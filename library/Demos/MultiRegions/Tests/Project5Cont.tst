<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Projection h-type convergence with 5 modes (continuous)</description>
    <executable>Project</executable>
    <parameters>Project_blank.xml -n 5</parameters>
    <files>
        <file description="Blank Session File">Project_blank.xml</file>
    </files>

    <metrics>
        <metric type="REGEX" id="1">
            <regex>
                Gradient: ([+-]?\d+(?:\.\d*)?)
            </regex>
            <matches>
                <match>
                    <field tolerance="1e-3" id="1">-4.9590018919</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>



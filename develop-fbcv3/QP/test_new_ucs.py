import cbmpy

# load model and old user constraints
ucmod = cbmpy.readSBML3FBC('test_user_def_rc.xml')
ucmod.importUserConstraints('test_user_def_rc.xml.user_constraints.json')
print(ucmod.user_constraints)

ucmod.getActiveObjective().getFluxObjectives()[0].setType('quadratic')

for r in ucmod.reactions[:3]:
    r.setNotes('<p>Finally <strong>notes</strong> are free to be notes!</p>')

for s in ucmod.species:
    if s.is_boundary:
        s.setAnnotation('secret', s.getId())

for c in ucmod.compartments:
    c.setAnnotation('secret', c.getId())

for p in ucmod.parameters:
    p.setAnnotation('secret', p.getId())

for o in ucmod.objectives:
    o.setAnnotation('secret', o.getId())
    for fo in o.flux_objectives:
        fo.setAnnotation('secret', fo.getId())

# convert to new FBCv3 structures
ucmod.convertUserConstraintsToUserDefinedConstraints()
print(ucmod.user_constraints, ucmod.user_defined_constraints)
print(list(ucmod.__global_id__.keys()))


for uc in ucmod.user_defined_constraints:
    uc.setAnnotation('secret', uc.getId())
    for var in uc.constraint_components:
        var.setAnnotation('secret', var.getId())


notes = """
        <html:p>
          <html:span style="font-family: Courier New,Courier,monospace;">
            <span xmlns="http://www.w3.org/1999/xhtml">
              <p>
                <span style="font-family: Courier New,Courier,monospace;"/>
              </p>
              <p>
                <br/>
                Core MEMESA model extended and updated as an FBC version 3 test model.<br/>
                Build with <a href="https://systemsbioinformatics.github.io/cbmpy/">CBMPy</a> 0.8.7
              </p>
            </span>
          </html:span>
        </html:p>

"""
ucmod.setNotes(notes)

# write my first SBML3 FBCv3 file

cbmpy.CBWrite.writeSBML3FBCV3(ucmod, 'the_third.fbc3.xml', compress_bounds=True)

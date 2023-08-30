import cbmpy

# load model and old user constraints
ucmod = cbmpy.readSBML3FBC('test_user_def_rc.xml')
ucmod.importUserConstraints('test_user_def_rc.xml.user_constraints.json')
print(ucmod.user_constraints)


for r in ucmod.reactions:
    r.setNotes('<p>Finally <strong>notes</strong> are free to be notes!</p>')

for s in ucmod.species:
    if s.is_boundary:
        s.setAnnotation('IsBoundaryOf', ucmod.getId())


# convert to new FBCv3 structures
ucmod.convertUserConstraintsToUserDefinedConstraints()
print(ucmod.user_constraints, ucmod.user_defined_constraints)
print(list(ucmod.__global_id__.keys()))

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

cbmpy.CBWrite.writeSBML3FBCV3(ucmod, 'the_third.v3.xml', compress_bounds=True)

import cbmpy

# load model and old user constraints
ucmod = cbmpy.readSBML3FBC('test_user_def_rc.xml')
ucmod.importUserConstraints('test_user_def_rc.xml.user_constraints.json')
#print(ucmod.user_constraints)

for r in ucmod.reactions[:3]:
    r.setNotes('<p>Finally <strong>notes</strong> are free to be notes!</p>')

for s in ucmod.species:
    if s.is_boundary:
        s.setAnnotation('secret', s.getId())

for c in ucmod.compartments:
    c.setAnnotation('secret', c.getId())

# ucmod.getActiveObjective().getFluxObjectives()[0].setType('quadratic')
ucmod.getActiveObjective().setAnnotation('quadratic_objective', '0.5*R01*R01,0.5*R24*R25')
for fo in ucmod.getActiveObjective().flux_objectives:
    fo.setAnnotation('secret', fo.getId())

print('\nFBA with OLD constraints:', cbmpy.doFBA(ucmod))
old_cons =  (ucmod.getReaction('R14').getValue(), ucmod.getReaction('R15').getValue(), ucmod.getReaction('R16').getValue())
print(old_cons)
print('\n')

# convert to new FBCv3 structures
ucmod.convertUserConstraintsToUserDefinedConstraints()
print('usercons', ucmod.user_constraints, ucmod.user_defined_constraints)
print('global_id', list(ucmod.__global_id__.keys()))

print('\nFBA with NEW constraints:', cbmpy.doFBA(ucmod))
new_cons = (ucmod.getReaction('R14').getValue(), ucmod.getReaction('R15').getValue(), ucmod.getReaction('R16').getValue())
print(new_cons)

print(old_cons[0] - new_cons[0] == 0)
print(old_cons[1] - new_cons[1] == 0)
print(old_cons[2] - new_cons[2] == 0)



for uc in ucmod.user_defined_constraints:
    uc.setAnnotation('secret', uc.getId())
    for var in uc.constraint_components:
        var.setAnnotation('secret', var.getId())


P1 = ucmod.createParameter('mypar1', 1)
P1.setAnnotation('secret', 'mypar1')
P2 = ucmod.createParameter('mypar2', 2)
P2.setAnnotation('secret', 'mypar2')
ucmod.addParameter(P1)
ucmod.addParameter(P2)



for p in ucmod.parameters:
    p.setAnnotation('secret', p.getId())


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

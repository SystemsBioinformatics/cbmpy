"""
MODEL DEVELOPMENT AREA

This module contains multiple environment development code for the
PySCeS Constraint Based Model package.
Brett G. Olivier (bgoli@users.sourceforge.net)

(c) Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam 2010.
All rights reserved.
"""

def Define_model_template():
    model_name = 'The model name'
    
    Reactions = {'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''}
                }
                
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True, 'SUBSYSTEM' : ''},
                'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'C1'}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'obj1' : {'id' : 'obj1', 'flux' : 'R01', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    print '\nModel name:', model_name
    print '\nReactions:'
    print Reactions
    print '\nSpecies:'
    print Species
    print '\nBounds:'
    print Bounds
    print '\nObjective function:'
    print Objective_function
    
    return None, None, None, None, None
    
def Define_core_model_1():
    """\nOriginal core model\n"""
    
    model_name = 'core_model_1'
    
    Reactions ={'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''},
                'R02' : {'id' : 'R02', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'C1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'C'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R05' : {'id' : 'R05', 'reversible' : False, 'reagents' : [(-1, 'B'), (1, 'D')], 'SUBSYSTEM' : ''},
                'R06' : {'id' : 'R06', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'E1')], 'SUBSYSTEM' : 'C2'},
                'R07' : {'id' : 'R07', 'reversible' : False, 'reagents' : [(-1, 'E1'), (1, 'E2')], 'SUBSYSTEM' : 'C2'},
                'R08' : {'id' : 'R08', 'reversible' : False, 'reagents' : [(-1, 'E2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R09' : {'id' : 'R09', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'F1')], 'SUBSYSTEM' : 'C2'},
                'R10' : {'id' : 'R10', 'reversible' : False, 'reagents' : [(-1, 'F1'), (1, 'F2')], 'SUBSYSTEM' : 'C2'},
                'R11' : {'id' : 'R11', 'reversible' : False, 'reagents' : [(-1, 'F2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R12' : {'id' : 'R12', 'reversible' : False, 'reagents' : [(-1, 'G'), (1, 'H')], 'SUBSYSTEM' : ''},
                'R13' : {'id' : 'R13', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'I1')], 'SUBSYSTEM' : 'C3'},
                'R14' : {'id' : 'R14', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R15' : {'id' : 'R15', 'reversible' : False, 'reagents' : [(-1, 'I2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R16' : {'id' : 'R16', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'J1')], 'SUBSYSTEM' : 'C3'},
                'R17' : {'id' : 'R17', 'reversible' : False, 'reagents' : [(-1, 'J1'), (1, 'J2')], 'SUBSYSTEM' : 'C3'},
                'R18' : {'id' : 'R18', 'reversible' : False, 'reagents' : [(-1, 'J2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R19' : {'id' : 'R19', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'K1')], 'SUBSYSTEM' : 'C3'},
                'R20' : {'id' : 'R20', 'reversible' : True, 'reagents' : [(-1, 'K1'), (1, 'K2')], 'SUBSYSTEM' : 'C3'},
                'R21' : {'id' : 'R21', 'reversible' : True, 'reagents' : [(-1, 'K2'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R22' : {'id' : 'R22', 'reversible' : False, 'reagents' : [(-1, 'L'), (1, 'M')], 'SUBSYSTEM' : ''},
                'R23' : {'id' : 'R23', 'reversible' : True, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R24' : {'id' : 'R24', 'reversible' : False, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R25' : {'id' : 'R25', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : ''}
               }
               
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True, 'SUBSYSTEM' : ''},
                'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'B' : {'id' : 'B', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'C' : {'id' : 'C', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'D' : {'id' : 'D', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E1' : {'id' : 'E1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E2' : {'id' : 'E2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F1' : {'id' : 'F1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F2' : {'id' : 'F2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'G' : {'id' : 'G', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'H' : {'id' : 'H', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I1' : {'id' : 'I1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I2' : {'id' : 'I2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J1' : {'id' : 'J1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J2' : {'id' : 'J2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K1' : {'id' : 'K1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K2' : {'id' : 'K2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'L' : {'id' : 'L', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'M' : {'id' : 'M', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'X1' : {'id' : 'X1', 'boundary' : True, 'SUBSYSTEM' : ''}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'objMaxJ25' : {'id' : 'objMaxJ25', 'flux' : 'R25', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_core_memesa_model():
    """\nOriginal core model + inefficient branch\n"""
    
    model_name = 'core_model_1b'
    
    Reactions ={'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''},
                'R02' : {'id' : 'R02', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'C1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'C'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R05' : {'id' : 'R05', 'reversible' : False, 'reagents' : [(-1, 'B'), (1, 'D')], 'SUBSYSTEM' : ''},
                'R06' : {'id' : 'R06', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'E1')], 'SUBSYSTEM' : 'C2'},
                'R07' : {'id' : 'R07', 'reversible' : False, 'reagents' : [(-1, 'E1'), (1, 'E2')], 'SUBSYSTEM' : 'C2'},
                'R08' : {'id' : 'R08', 'reversible' : False, 'reagents' : [(-1, 'E2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R09' : {'id' : 'R09', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'F1')], 'SUBSYSTEM' : 'C2'},
                'R10' : {'id' : 'R10', 'reversible' : False, 'reagents' : [(-1, 'F1'), (1, 'F2')], 'SUBSYSTEM' : 'C2'},
                'R11' : {'id' : 'R11', 'reversible' : False, 'reagents' : [(-1, 'F2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R12' : {'id' : 'R12', 'reversible' : False, 'reagents' : [(-1, 'G'), (1, 'H')], 'SUBSYSTEM' : ''},
                'R13' : {'id' : 'R13', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'I1')], 'SUBSYSTEM' : 'C3'},
                'R14' : {'id' : 'R14', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R15' : {'id' : 'R15', 'reversible' : False, 'reagents' : [(-1, 'I2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R16' : {'id' : 'R16', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'J1')], 'SUBSYSTEM' : 'C3'},
                'R17' : {'id' : 'R17', 'reversible' : False, 'reagents' : [(-1, 'J1'), (1, 'J2')], 'SUBSYSTEM' : 'C3'},
                'R18' : {'id' : 'R18', 'reversible' : False, 'reagents' : [(-1, 'J2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R19' : {'id' : 'R19', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'K1')], 'SUBSYSTEM' : 'C3'},
                'R20' : {'id' : 'R20', 'reversible' : True, 'reagents' : [(-1, 'K1'), (1, 'K2')], 'SUBSYSTEM' : 'C3'},
                'R21' : {'id' : 'R21', 'reversible' : True, 'reagents' : [(-1, 'K2'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R22' : {'id' : 'R22', 'reversible' : False, 'reagents' : [(-1, 'L'), (1, 'M')], 'SUBSYSTEM' : ''},
                'R23' : {'id' : 'R23', 'reversible' : True, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R24' : {'id' : 'R24', 'reversible' : False, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R25' : {'id' : 'R25', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : ''},
                'R26' : {'id' : 'R26', 'reversible' : False, 'reagents' : [(-1, 'A'), (0.5, 'N'), (0.5, 'X3')], 'SUBSYSTEM' : ''}
               }
               
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True, 'SUBSYSTEM' : ''},
                'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'B' : {'id' : 'B', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'C' : {'id' : 'C', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'D' : {'id' : 'D', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E1' : {'id' : 'E1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E2' : {'id' : 'E2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F1' : {'id' : 'F1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F2' : {'id' : 'F2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'G' : {'id' : 'G', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'H' : {'id' : 'H', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I1' : {'id' : 'I1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I2' : {'id' : 'I2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J1' : {'id' : 'J1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J2' : {'id' : 'J2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K1' : {'id' : 'K1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K2' : {'id' : 'K2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'L' : {'id' : 'L', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'M' : {'id' : 'M', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'X1' : {'id' : 'X1', 'boundary' : True, 'SUBSYSTEM' : ''},
                'X3' : {'id' : 'X3', 'boundary' : True, 'SUBSYSTEM' : ''}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'objMaxJ25' : {'id' : 'objMaxJ25', 'flux' : 'R25', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function

    
def Define_core_model_2a():
    """\nCore model with all reactions reversible\n"""
    
    model_name = 'core_model_2a'
    
    Reactions ={'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''},
                'R02' : {'id' : 'R02', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'C1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'C'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R05' : {'id' : 'R05', 'reversible' : True, 'reagents' : [(-1, 'B'), (1, 'D')], 'SUBSYSTEM' : ''},
                'R06' : {'id' : 'R06', 'reversible' : True, 'reagents' : [(-1, 'D'), (1, 'E1')], 'SUBSYSTEM' : 'C2'},
                'R07' : {'id' : 'R07', 'reversible' : True, 'reagents' : [(-1, 'E1'), (1, 'E2')], 'SUBSYSTEM' : 'C2'},
                'R08' : {'id' : 'R08', 'reversible' : True, 'reagents' : [(-1, 'E2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R09' : {'id' : 'R09', 'reversible' : True, 'reagents' : [(-1, 'D'), (1, 'F1')], 'SUBSYSTEM' : 'C2'},
                'R10' : {'id' : 'R10', 'reversible' : True, 'reagents' : [(-1, 'F1'), (1, 'F2')], 'SUBSYSTEM' : 'C2'},
                'R11' : {'id' : 'R11', 'reversible' : True, 'reagents' : [(-1, 'F2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R12' : {'id' : 'R12', 'reversible' : True, 'reagents' : [(-1, 'G'), (1, 'H')], 'SUBSYSTEM' : ''},
                'R13' : {'id' : 'R13', 'reversible' : True, 'reagents' : [(-1, 'H'), (1, 'I1')], 'SUBSYSTEM' : 'C3'},
                'R14' : {'id' : 'R14', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R15' : {'id' : 'R15', 'reversible' : True, 'reagents' : [(-1, 'I2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R16' : {'id' : 'R16', 'reversible' : True, 'reagents' : [(-1, 'H'), (1, 'J1')], 'SUBSYSTEM' : 'C3'},
                'R17' : {'id' : 'R17', 'reversible' : True, 'reagents' : [(-1, 'J1'), (1, 'J2')], 'SUBSYSTEM' : 'C3'},
                'R18' : {'id' : 'R18', 'reversible' : True, 'reagents' : [(-1, 'J2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R19' : {'id' : 'R19', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'K1')], 'SUBSYSTEM' : 'C3'},
                'R20' : {'id' : 'R20', 'reversible' : True, 'reagents' : [(-1, 'K1'), (1, 'K2')], 'SUBSYSTEM' : 'C3'},
                'R21' : {'id' : 'R21', 'reversible' : True, 'reagents' : [(-1, 'K2'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R22' : {'id' : 'R22', 'reversible' : True, 'reagents' : [(-1, 'L'), (1, 'M')], 'SUBSYSTEM' : ''},
                'R23' : {'id' : 'R23', 'reversible' : True, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R24' : {'id' : 'R24', 'reversible' : True, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R25' : {'id' : 'R25', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : ''}
               }
               
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True, 'SUBSYSTEM' : ''},
                'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'B' : {'id' : 'B', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'C' : {'id' : 'C', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'D' : {'id' : 'D', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E1' : {'id' : 'E1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E2' : {'id' : 'E2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F1' : {'id' : 'F1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F2' : {'id' : 'F2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'G' : {'id' : 'G', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'H' : {'id' : 'H', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I1' : {'id' : 'I1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I2' : {'id' : 'I2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J1' : {'id' : 'J1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J2' : {'id' : 'J2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K1' : {'id' : 'K1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K2' : {'id' : 'K2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'L' : {'id' : 'L', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'M' : {'id' : 'M', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'X1' : {'id' : 'X1', 'boundary' : True, 'SUBSYSTEM' : ''}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'obj1' : {'id' : 'obj1', 'flux' : 'R25', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_core_model_3():
    """\nCore model with added 'cofactor' cycle\n"""
    
    model_name = 'core_model_3'
    
    Reactions ={'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''},
                'R02' : {'id' : 'R02', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'C1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'C'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R05' : {'id' : 'R05', 'reversible' : False, 'reagents' : [(-0.5, 'B'), (1, 'D')], 'SUBSYSTEM' : ''},
                'R06' : {'id' : 'R06', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'E1')], 'SUBSYSTEM' : 'C2'},
                'R07' : {'id' : 'R07', 'reversible' : False, 'reagents' : [(-1, 'E1'), (1, 'E2')], 'SUBSYSTEM' : 'C2'},
                'R08' : {'id' : 'R08', 'reversible' : False, 'reagents' : [(-1, 'E2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R09' : {'id' : 'R09', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'F1')], 'SUBSYSTEM' : 'C2'},
                'R10' : {'id' : 'R10', 'reversible' : False, 'reagents' : [(-1, 'F1'), (1, 'F2')], 'SUBSYSTEM' : 'C2'},
                'R11' : {'id' : 'R11', 'reversible' : False, 'reagents' : [(-1, 'F2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
                'R12' : {'id' : 'R12', 'reversible' : False, 'reagents' : [(-1, 'G'), (1, 'H')], 'SUBSYSTEM' : ''},
                'R13' : {'id' : 'R13', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'I1')], 'SUBSYSTEM' : 'C3'},
                'R14' : {'id' : 'R14', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                'R15' : {'id' : 'R15', 'reversible' : False, 'reagents' : [(-1, 'I2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R16' : {'id' : 'R16', 'reversible' : False, 'reagents' : [(-1, 'H'), (1, 'J1')], 'SUBSYSTEM' : 'C3'},
                'R17' : {'id' : 'R17', 'reversible' : False, 'reagents' : [(-1, 'J1'), (1, 'J2')], 'SUBSYSTEM' : 'C3'},
                'R18' : {'id' : 'R18', 'reversible' : False, 'reagents' : [(-1, 'J2'), (1, 'L')], 'SUBSYSTEM' : 'C3'},
                'R19' : {'id' : 'R19', 'reversible' : True, 'reagents' : [(-1, 'I1'), (1, 'K1')], 'SUBSYSTEM' : 'C3'},
                'R20' : {'id' : 'R20', 'reversible' : True, 'reagents' : [(-1, 'K1'), (1, 'K2')], 'SUBSYSTEM' : 'C3'},
                'R21' : {'id' : 'R21', 'reversible' : True, 'reagents' : [(-1, 'K2'), (1, 'I2')], 'SUBSYSTEM' : 'C3'},
                ##  'R22' : {'id' : 'R22', 'reversible' : False, 'reagents' : [(-1, 'L'), (1, 'M')], 'SUBSYSTEM' : ''},
                'R23' : {'id' : 'R23', 'reversible' : True, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R24' : {'id' : 'R24', 'reversible' : False, 'reagents' : [(-1, 'M'), (1, 'N')], 'SUBSYSTEM' : 'C4'},
                'R25' : {'id' : 'R25', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : ''},
                'R26' : {'id' : 'R26', 'reversible' : False, 'reagents' : [(-0.5, 'B'), (-1, 'X2'), (1, 'P')], 'SUBSYSTEM' : 'E1'},
                'R27' : {'id' : 'R27', 'reversible' : True, 'reagents' : [(-1, 'P'), (1, 'R')], 'SUBSYSTEM' : 'E1'},
                'R28' : {'id' : 'R28', 'reversible' : False, 'reagents' : [(-1, 'L'), (-1, 'R'), (1, 'M'), (1, 'X3')], 'SUBSYSTEM' : 'E1'}
               }
               
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True, 'SUBSYSTEM' : ''},
                'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'B' : {'id' : 'B', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'C' : {'id' : 'C', 'boundary' : False, 'SUBSYSTEM' : 'C1'},
                'D' : {'id' : 'D', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E1' : {'id' : 'E1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'E2' : {'id' : 'E2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F1' : {'id' : 'F1', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'F2' : {'id' : 'F2', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'G' : {'id' : 'G', 'boundary' : False, 'SUBSYSTEM' : 'C2'},
                'H' : {'id' : 'H', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I1' : {'id' : 'I1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'I2' : {'id' : 'I2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J1' : {'id' : 'J1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'J2' : {'id' : 'J2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K1' : {'id' : 'K1', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'K2' : {'id' : 'K2', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'L' : {'id' : 'L', 'boundary' : False, 'SUBSYSTEM' : 'C3'},
                'M' : {'id' : 'M', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'C4'},
                'X1' : {'id' : 'X1', 'boundary' : True, 'SUBSYSTEM' : ''},
                'P' : {'id' : 'P', 'boundary' : False, 'SUBSYSTEM' : 'E1'},
                'R' : {'id' : 'R', 'boundary' : False, 'SUBSYSTEM' : 'E1'},
                'X2' : {'id' : 'X2', 'boundary' : True, 'SUBSYSTEM' : 'E1'},
                'X3' : {'id' : 'X3', 'boundary' : True, 'SUBSYSTEM' : 'E1'}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'obj1' : {'id' : 'obj1', 'flux' : 'R25', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_milp_model_1():
    """\nOriginal MILP model\n"""
    
    model_name = 'core_model_1'
    
    Reactions ={'RA' : {'id' : 'RA', 'reversible' : True, 'reagents' : [(1, 'A')], 'SUBSYSTEM' : 'b1'},
                'RB' : {'id' : 'RB', 'reversible' : True, 'reagents' : [(1, 'B')], 'SUBSYSTEM' : 'b2'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'b1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'B'), (1, 'C')], 'SUBSYSTEM' : 'b2'},
                'R05' : {'id' : 'R05', 'reversible' : False, 'reagents' : [(-1, 'C')], 'SUBSYSTEM' : 'b3'}
               }
               
    Species = { 'A' : {'id' : 'A', 'boundary' : False, 'SUBSYSTEM' : 'b1'},
                'B' : {'id' : 'B', 'boundary' : False, 'SUBSYSTEM' : 'b2'},
                'C' : {'id' : 'C', 'boundary' : False, 'SUBSYSTEM' : 'b3'}
              }
              
    Bounds = {'RA' : {'lower' : 2, 'upper' : 10},
              'RB' : {'lower' : -10, 'upper' : 0},
              'R03' : {'lower' : -10, 'upper' : 10},
              'R04' : {'lower' : -10, 'upper' : 10},
              'R05' : {'lower' : 0, 'upper' : 0}
              }
    
    Objective_function = {'objMaxR05' : {'id' : 'objMaxR05', 'flux' : 'R05', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_flp_model_1():
    """\nThe fractional programming model\n"""
    
    model_name = 'frac_model_1'
    
    Reactions ={
                'R0' : {'id' : 'R0', 'reversible' : False, 'reagents' : [(-1, 'N0'), (1, 'N')], 'SUBSYSTEM' : 'G'},
                'R1' : {'id' : 'R1', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : 'C'},
                'R5' : {'id' : 'R5', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X3')], 'SUBSYSTEM' : 'A'},
                'R4' : {'id' : 'R4', 'reversible' : False, 'reagents' : [(-1, 'X2'), (-2, 'X4'), (-2, 'ATP'), (1, 'NX')], 'SUBSYSTEM' : 'G'},
                'R2' : {'id' : 'R2', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (4, 'ATP')], 'SUBSYSTEM' : 'C'},
                'R3' : {'id' : 'R3', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (2, 'ATP')], 'SUBSYSTEM' : 'C'},
                'R6' : {'id' : 'R6', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4'), (-2, 'ATP')], 'SUBSYSTEM' : 'A'},
                'R7' : {'id' : 'R7', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4')], 'SUBSYSTEM' : 'A'},
               }
               
    Species = {
                'N0' : {'id' : 'N0', 'boundary' : True, 'SUBSYSTEM' : 'A'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X1' : {'id' : 'X1', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X2' : {'id' : 'X2', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X3' : {'id' : 'X3', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'X4' : {'id' : 'X4', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'ATP' : {'id' : 'ATP', 'boundary' : False, 'SUBSYSTEM' : 'G'},
                'NX' : {'id' : 'NX', 'boundary' : True, 'SUBSYSTEM' : 'G'},
              }
              
    Bounds = {
                'R0' : {'lower' : 0, 'upper' : 10.0},

              }

    
    Objective_function = {'objMaxR4' : {'id' : 'objMaxR4', 'flux' : 'R4', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
    
def Define_flp_model_1b():
    """\nThe fractional programming model with ATPase\n"""
    
    model_name = 'frac_model_1b'
    
    Reactions ={
                'R0' : {'id' : 'R0', 'reversible' : False, 'reagents' : [(-1, 'N0'), (1, 'N')], 'SUBSYSTEM' : 'G'},
                'R1' : {'id' : 'R1', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : 'C'},
                'R5' : {'id' : 'R5', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X3')], 'SUBSYSTEM' : 'A'},
                'R4' : {'id' : 'R4', 'reversible' : False, 'reagents' : [(-1, 'X2'), (-2, 'X4'), (-2, 'ATP'), (1, 'NX')], 'SUBSYSTEM' : 'G'},
                'R2' : {'id' : 'R2', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (4, 'ATP')], 'SUBSYSTEM' : 'C'},
                'R3' : {'id' : 'R3', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (2, 'ATP')], 'SUBSYSTEM' : 'C'},
                'R6' : {'id' : 'R6', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4'), (-2, 'ATP')], 'SUBSYSTEM' : 'A'},
                'R7' : {'id' : 'R7', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4')], 'SUBSYSTEM' : 'A'},
                'RATPase' : {'id' : 'RATPase', 'reversible' : False, 'reagents' : [(-1, 'ATP'), (1, 'P')], 'SUBSYSTEM' : 'G'},
                'RATPsyn' : {'id' : 'RATPsyn', 'reversible' : False, 'reagents' : [(1, 'ATP'), (-1, 'P')], 'SUBSYSTEM' : 'G'}
               }
               
    Species = {
                'N0' : {'id' : 'N0', 'boundary' : True, 'SUBSYSTEM' : 'A'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X1' : {'id' : 'X1', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X2' : {'id' : 'X2', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X3' : {'id' : 'X3', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'X4' : {'id' : 'X4', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'ATP' : {'id' : 'ATP', 'boundary' : False, 'SUBSYSTEM' : 'G'},
                'NX' : {'id' : 'NX', 'boundary' : True, 'SUBSYSTEM' : 'G'},
                'P' : {'id' : 'P', 'boundary' : False, 'SUBSYSTEM' : 'G'},
              }
              
    Bounds = {
                'R0' : {'lower' : 0, 'upper' : 10.0},
              }

    
    Objective_function = {'objMaxR4' : {'id' : 'objMaxR4', 'flux' : 'R4', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_flp_model_2():
    """\nThe fractional programming model with ATPase and ADP\n"""
    
    model_name = 'frac_model_2'
    
    Reactions ={
                'R0' : {'id' : 'R0', 'reversible' : False, 'reagents' : [(-1, 'N0'), (1, 'N')], 'SUBSYSTEM' : 'G'},
                'R1' : {'id' : 'R1', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X1')], 'SUBSYSTEM' : 'C'},
                'R5' : {'id' : 'R5', 'reversible' : False, 'reagents' : [(-1, 'N'), (1, 'X3')], 'SUBSYSTEM' : 'A'},
                'R4' : {'id' : 'R4', 'reversible' : False, 'reagents' : [(-1, 'X2'), (-2, 'X4'), (-2, 'ATP'), (-2, 'ADP'), (1, 'NX')], 'SUBSYSTEM' : 'G'},
                'R2' : {'id' : 'R2', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (4, 'ATP'), (-4, 'ADP')], 'SUBSYSTEM' : 'C'},
                'R3' : {'id' : 'R3', 'reversible' : False, 'reagents' : [(-1, 'X1'), (1, 'X2'), (2, 'ATP'), (-2, 'ADP')], 'SUBSYSTEM' : 'C'},
                'R6' : {'id' : 'R6', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4'), (-2, 'ATP'), (2, 'ADP')], 'SUBSYSTEM' : 'A'},
                'R7' : {'id' : 'R7', 'reversible' : False, 'reagents' : [(-1, 'X3'), (1, 'X4')], 'SUBSYSTEM' : 'A'},
                'RATPase' : {'id' : 'RATPase', 'reversible' : False, 'reagents' : [(-1, 'ATP'), (1, 'ADP'), (1, 'P')], 'SUBSYSTEM' : 'G'},
                'RATPsyn' : {'id' : 'RATPsyn', 'reversible' : False, 'reagents' : [(1, 'ATP'), (-1, 'ADP'), (-1, 'P')], 'SUBSYSTEM' : 'G'}
               }
               
    Species = {
                'N0' : {'id' : 'N0', 'boundary' : True, 'SUBSYSTEM' : 'A'},
                'N' : {'id' : 'N', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X1' : {'id' : 'X1', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X2' : {'id' : 'X2', 'boundary' : False, 'SUBSYSTEM' : 'A'},
                'X3' : {'id' : 'X3', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'X4' : {'id' : 'X4', 'boundary' : False, 'SUBSYSTEM' : 'C'},
                'ATP' : {'id' : 'ATP', 'boundary' : False, 'SUBSYSTEM' : 'G'},
                'NX' : {'id' : 'NX', 'boundary' : True, 'SUBSYSTEM' : 'G'},
                'P' : {'id' : 'P', 'boundary' : False, 'SUBSYSTEM' : 'G'},
              }
              
    Bounds = {
                'R0' : {'lower' : 0, 'upper' : 10.0},
              }
    
    Objective_function = {'objMaxR4' : {'id' : 'objMaxR4', 'flux' : 'R4', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function
    
def Define_gpr_test_model():
    """\nModel to test various GPR encodings\n"""
    
    model_name = 'test_gpr'
    
    Reactions ={'R01' : {'id' : 'R01', 'reversible' : False, 'reagents' : [(-1, 'X0'), (1, 'A')], 'SUBSYSTEM' : ''},
                'R02' : {'id' : 'R02', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R03' : {'id' : 'R03', 'reversible' : True, 'reagents' : [(-1, 'A'), (1, 'C')], 'SUBSYSTEM' : 'C1'},
                'R04' : {'id' : 'R04', 'reversible' : True, 'reagents' : [(-1, 'C'), (1, 'B')], 'SUBSYSTEM' : 'C1'},
                'R05' : {'id' : 'R05', 'reversible' : False, 'reagents' : [(-1, 'B'), (1, 'D')], 'SUBSYSTEM' : ''},
                'R06' : {'id' : 'R06', 'reversible' : False, 'reagents' : [(-1, 'D'), (1, 'E1')], 'SUBSYSTEM' : 'C2'},
                'R07' : {'id' : 'R07', 'reversible' : False, 'reagents' : [(-1, 'E1'), (1, 'E2')], 'SUBSYSTEM' : 'C2'},
                'R08' : {'id' : 'R08', 'reversible' : False, 'reagents' : [(-1, 'E2'), (1, 'G')], 'SUBSYSTEM' : 'C2'},
               }
               
    Species = { 'X0' : {'id' : 'X0', 'boundary' : True},
                'A' : {'id' : 'A', 'boundary' : False},
                'B' : {'id' : 'B', 'boundary' : False},
                'C' : {'id' : 'C', 'boundary' : False},
                'D' : {'id' : 'D', 'boundary' : False},
                'E1' : {'id' : 'E1', 'boundary' : False},
                'E2' : {'id' : 'E2', 'boundary' : False},
                'G' : {'id' : 'G', 'boundary' : True}
              }
              
    Bounds = {'R01' : {'lower' : 0, 'upper' : 1}}
    
    Objective_function = {'obj1' : {'id' : 'obj1', 'flux' : 'R08', 'coefficient' : 1, 'sense' : 'Maximize', 'active' : True}}
    
    return model_name, Reactions, Species, Bounds, Objective_function

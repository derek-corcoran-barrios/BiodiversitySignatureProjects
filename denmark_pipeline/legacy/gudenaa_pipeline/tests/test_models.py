"""Small, deliberately uneven subareas test model policy and checkpoint syntax.
Requires amplpy plus AMPL base/HiGHS modules. Run from the project root.
These tests do not substitute for a full Gudenå spatial run.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from amplpy import AMPL

ROOT = Path(__file__).resolve().parents[1]
H = ['ForestDryPoor', 'ForestDryRich', 'ForestWetPoor', 'ForestWetRich',
     'OpenDryPoor', 'OpenDryRich', 'OpenWetPoor', 'OpenWetRich']
K = ['Global', 'Subarea1', 'Subarea2', 'Subarea3', 'Subarea4']
CELLS = list(range(1, 5))
PAIRS = [(h, c) for h in H for c in CELLS if c != 1 or h.startswith('Forest')]

def load(ampl, model):
    ampl.read(str(ROOT / 'models' / model))
    for name, value in [('Cells', CELLS), ('Landuses', H), ('Subareas', CELLS),
                        ('AllowedActions', PAIRS), ('E', [(1, 2), (2, 3), (3, 4)]),
                        ('SameHabitatEdges', [])]:
        ampl.get_set(name).set_values(value)
    params = {'BioDiversity': {p: H.index(p[0])+1 for p in PAIRS},
              'RestorableAreaHa': {c: 100 for c in CELLS},
              'RestorableAreaHaBySubarea': {(c,c): 100 for c in CELLS},
              'SubareaAreaHa': dict(zip(CELLS, [2000,8000,2000,2000])),
              'CurrentNatureAreaHa': {c: 950 for c in CELLS},
              'FinalNatureTargetHa': dict(zip(CELLS,[600,2400,600,600])),
              'RestorationRequirementHa': {c: 0 for c in CELLS},
              'BiodiversityNormalizer': 32, 'ContiguityNormalizer': 3,
              'ExistingNatureContactNormalizer': 1, 'SpatialContiguityBonus': 1.2,
              'NewToNewAdjacencyWeight': .5, 'ExistingToNewAdjacencyWeight': .5,
              'GlobalRestorationToleranceHa': 0}
    if model.startswith('network'):
        for name, value in [('ForestLanduses',H[:4]), ('RequiredNetworks',K),
            ('GlobalNetworks',['Global']), ('ExistingPatches',[1]),
            ('CellPatchEdges',[(c,1) for c in CELLS]),
            ('EligibleNetworkCells',[(k,c) for k in K for c in CELLS]),
            ('EligibleNetworkPatches',[(k,1) for k in K]),
            ('NetworkCellArcs',[]), ('NetworkCellPatchEdges',[(k,c,1) for k in K for c in CELLS])]:
            ampl.get_set(name).set_values(value)
        params.update({'GlobalPolygonAreaHa':14000, 'GlobalCurrentNatureAreaHa':3800,
            'GlobalFinalNatureTargetHa':4200, 'GlobalRestorationRequirementHa':400,
            'ExistingPatchAreaHa':{1:4800}, 'NetworkMinimumHa':dict(zip(K,[5000,1000,1000,1000,1000])),
            'NetworkCellAreaHa':{(k,c):100 for k in K for c in CELLS if k=='Global' or k==f'Subarea{c}'},
            'NetworkPatchAreaHa':{(k,1):4800 if k=='Global' else 950 for k in K},
            'ComponentMaxAreaCellHa':{c:5200 for c in CELLS},
            'ComponentMaxAreaPatchHa':{1:5200}, 'IsProductionForestCell':{1:1},
            'FootprintChangeLimitHa':400})
    else:
        params['RestorationToleranceHa']={c:0 for c in CELLS}
    for name, value in params.items():
        ampl.get_parameter(name).set_values(value) if isinstance(value,dict) else ampl.get_parameter(name).set(value)
    ampl.option['solver']='highs'
    ampl.option['highs_options']='outlev=0 mip:bestbound=1'

def solve(ampl):
    ampl.solve()
    assert ampl.get_value('solve_result_num') < 100

def main():
    a=AMPL(); load(a,'warmup_global30.mod'); solve(a)
    assert a.get_value('sum {c in Cells} 100*Selected[c]') == 400
    # Subarea 2 ends with 13.125% nature, which is legal under the global policy.
    assert (950 + a.get_value('Selected[2]')*100)/8000 < .30
    assert a.get_value("LanduseDecision['ForestWetRich',1]") == 1
    a.close()
    a=AMPL(); load(a,'network_global30.mod')
    a.eval('drop WarmFootprintChange; drop RepairQuality;')
    solve(a)
    assert a.get_value('RelativeNetworkShortfall') < 1e-9
    a.eval('fix {k in RequiredNetworks} NetworkShortfallHa[k] := 0; drop RelativeNetworkShortfall; restore RepairQuality;')
    solve(a)
    a.eval('if max {i in 1.._ncons} max(0,-_con[i].slack) > 1e-5 then { exit 1; };')
    a.eval('if max {j in 1.._nvars} max(0,-_var[j].slack) > 1e-5 then { exit 1; };')
    assert a.get_value('sum {(k,p) in EligibleNetworkPatches} NetworkRoot[k,p]') == 5
    with TemporaryDirectory() as td:
        p=Path(td)/'checkpoint.run'
        a.eval(f'printf {{j in 1.._nvars}} "let %s := %.17g;\\n", _varname[j], _var[j] > "{p}";')
        original=a.get_value('RepairQuality')
        b=AMPL(); load(b,'network_global30.mod'); b.read(str(p))
        assert abs(b.get_value('RepairQuality')-original) < 1e-9
        b.close()
    # An inadequate restoration budget produces a strictly positive optimum,
    # even though every potential component has enough total area.
    a.eval('unfix {k in RequiredNetworks} NetworkShortfallHa[k]; drop RepairQuality; restore RelativeNetworkShortfall;')
    a.get_parameter('GlobalPolygonAreaHa').set(4100/.3)
    a.get_parameter('SubareaAreaHa').set_values({1:2000,2:4100/.3-6000,3:2000,4:2000})
    a.get_parameter('GlobalFinalNatureTargetHa').set(4100)
    a.get_parameter('GlobalRestorationRequirementHa').set(300)
    solve(a)
    assert a.get_value('RelativeNetworkShortfall') > 0
    a.close()
    print('PASS: global budget, uneven subareas, forest actions, all five certificates, checkpoint round-trip, positive shortfall.')

if __name__=='__main__': main()

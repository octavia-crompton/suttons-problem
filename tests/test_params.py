from sutton import Params

def test_params_to_dict_has_expected_keys():
    p = Params()
    d = p.to_dict()
    for k in ["x","z","nx","nz","ustar_f","ustar_c","Q_c","Q_f","Q_a"]:
        assert k in d

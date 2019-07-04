def color_lookup(ecv_name):
    color_dict = {
        "default": dict({"default": "binary"}),
        "example": dict({"Data": "RdYlGn",
                         "Sequential":"YlGn",
                         "Diverging":"BrBG"}),
        "bdalb" : dict({"Data": "YlGn_r",
                         "Sequential":"PuBu",
                         "Diverging":"BrBG"}),
        "bhalb" : dict({"Data": "YlGn_r",
                         "Sequential":"PuBu",
                         "Diverging":"BrBG"}),
        "albdiffbnd": dict(),
        "albdirbnd": dict(),
        "albisccp": dict({"Data": "RdYlGn",
                         "Sequential":"YlGn",
                         "Diverging":"BrBG"}),
        "burntArea": dict(),
        "cct": dict(),
        "cfc11": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "cfc12": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "chldiat": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "chldiaz": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}), 
        "chlmisc": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "chlpico": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "clt": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "cod": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "hfds": dict(),
        "hus": dict(),
        "huss": dict({"Data": "GnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG"}),
        "hur": dict(),
        "hurs": dict({"Data": "GnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG"}),
        "intdic": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "lai": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBg"}),
        "lwp": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "mrsos": dict(),
        "od550aer": dict(),
        "phnat": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "pr": dict({"Data": "GnBu",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "reffclwtop": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "rhs": dict({"Data": "GnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG"}),
        "rlds": dict(),
        "rldscs": dict(),
        "rlus": dict(),
        "rluscs": dict(),
        "rlut": dict(),
        "rsds": dict(),
        "rsdscs": dict(),
        "rsdt": dict(),
        "rsus": dict(),
        "rsuscs": dict(),
        "rsut": dict(),
        "sf6": dict(),
        "sfcWind": dict({"Data": "YlGnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG_r"}),
        "sftgif": dict(),
        "sic": dict({"Data": "Blues_r",
                     "Sequential":"Blues_r",
                     "Diverging":"RdBu_r"}),
        "sidivvel": dict(),
        "sidmassdyn": dict(),
        "sidmassth": dict(),
        "siextentn": dict(),
        "siextents": dict(),
        "sispeed": dict({"Data": "cubehelix",
                     "Sequential":"cubehelix",
                     "Diverging":"RdBu_r"}),
        "sit": dict({"Data": "gist_earth",
                     "Sequential":"gist_earth",
                     "Diverging":"RdBu_r"}),
        "siu": dict({"Data": "cubehelix",
                     "Sequential":"cubehelix",
                     "Diverging":"RdBu_r"}),
        "siv": dict({"Data": "cubehelix",
                     "Sequential":"cubehelix",
                     "Diverging":"RdBu_r"}),
        "sm": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "snc": dict(),
        "snd": dict({"Data": "jet", #bone
                    "Sequential": "BuPu",
                    "Diverging": "Spectral"}),
        "swe": dict(),
        "so": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "sos": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "spco2": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "ta": dict({"Data": "RdYlGn",
                    "Sequential":"YlGn",
                    "Diverging":"BrBG"}),
        "talknat": dict({"Data": "viridis",
                       "Sequential":"cool",
                       "Diverging":"bwr"}),
        "tas": dict({"Data": "YlOrRd",
                    "Sequential":"Wistia",
                    "Diverging":"RdBu_r"}),
        "tasmax": dict(),
        "tasmin": dict(),
        "tauuo": dict(),
        "tauvo": dict(),
        "to": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "tos": dict(),
        "toz": dict(),
        "tpf": dict(),
        "ts": dict(),
        "ua": dict(),
        "uas": dict({"Data": "YlGnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG_r"}),
        "uo": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "uos": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "va": dict(),
        "vas": dict({"Data": "YlGnBu",
                    "Sequential":"Wistia",
                    "Diverging":"BrBG_r"}),
        "vmrch4": dict({"Data": "RdYlGn",
                        "Sequential":"YlGn",
                        "Diverging":"BrBG"}),
        "vmrco2": dict({"Data": "RdYlGn",
                        "Sequential":"YlGn",
                        "Diverging":"BrBG"}),
        "vmrstrato3": dict(),
        "vmro3": dict(),
        "vo": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "vos": dict({"Data": "viridis",
                    "Sequential":"YlGn",
                    "Diverging":"coolwarm"}),
        "xco2": dict({"Data": "RdYlGn_r",
                        "Sequential":"YlGn",
                        "Diverging":"BrBG"}),
        "xch4": dict({"Data": "RdYlGn_r",
                        "Sequential":"YlGn",
                        "Diverging":"BrBG"}),
        "zos": dict({"Data": "viridis",
                    "Sequential":"GnBu",
                    "Diverging":"RdBu_r"}),
    }
    if ecv_name in color_dict.keys():
        if color_dict[ecv_name]:
            return color_dict[ecv_name]
        else:
            return color_dict["default"]
    else:
        return color_dict["default"]
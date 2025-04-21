

class ERV_module():
    def __init__(self, eff_sen = 0.8, eff_lat=0.7):
        self.eff_sen = eff_sen
        self.eff_lat = eff_lat

    def calc_air(self, ra_t, ra_w, oa_t, oa_w):
        t_sa = self.eff_sen * abs(oa_t - ra_t) + min([ra_t, oa_t])
        w_sa = self.eff_lat * abs(oa_w - ra_w) + min([ra_w, oa_w])
        return t_sa, w_sa
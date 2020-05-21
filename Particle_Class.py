class SpinParticle:
    def __init__(self, spin, bond, placement):
        self.spin      = spin
        self.bond      = bond
        self.placement = placement

    def bond_strenght(self):
        return self.bond

from ..echo_pro import EchoPro


class Visualize:
    """
    A Class that creates all plots that
    are requested by the user

    Parameters
    ----------
    echopro : EchoPro
        EchoPro object that contains all necessary parameters
    """

    def __init__(self,
                 echopro: EchoPro):

        if echopro.bootstrapping_performed:
            print("Plotting ...")
        else:
            raise RuntimeError("Bootstrapping must be performed first!")
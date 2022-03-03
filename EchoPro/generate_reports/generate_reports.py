from ..echo_pro import EchoPro


class GenerateReports:
    """
    A Class that writes requested variables to
    consolidated files.

    Parameters
    ----------
    echopro : EchoPro
        EchoPro object that contains all necessary parameters
    """

    def __init__(self,
                 echopro: EchoPro):

        if echopro.bootstrapping_performed:
            print("Generate all consolidated files")
        else:
            raise RuntimeError("Bootstrapping must be performed first!")

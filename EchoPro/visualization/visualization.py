from ..survey import Survey


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
                 survey: Survey):

        if survey.bootstrapping_performed:
            print("Plotting ...")
        else:
            raise RuntimeError("Bootstrapping must be performed first!")
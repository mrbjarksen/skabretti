from manimlib.imports import *

class FirstScene(Scene):
    def construct(self):
        eq = TexMobject("\\int_{-\\infty}^\\infty e^{-x^2} \\, dx \\, = \\, \\sqrt{\\pi}")
        eq.scale(2)
        self.play(Write(eq))
        self.wait(2)
        eq.to_edge(RIGHT)
        self.wait(1)
        eq.to_corner(DL)
        self.wait(1)
        eq.move_to(np.array([0, 3, 0]))
        self.wait(1)
        eq.move_to(np.array([0, 0, 0]))
        self.wait(3)

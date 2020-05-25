from manimlib.imports import *


class Kubbur(Rectangle):
    CONFIG = {
        "mass": 1,
        "velocity": 0,
        "acceleration": 0,
        "fill_opacity": 1.0,
        "fill_color": WHITE,
        "stroke_opacity": 0.0,
        "show_label": True,
        "label": "m",
        "label_scale": 1,
        "label_color": BLACK,
        "width_to_height_ratio": (1 + 5 ** 0.5) / 2,
        "scale_with_mass": True,
    }

    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        self.height = 1
        if self.scale_with_mass:
            self.height *= self.mass ** 0.5
        self.width = self.height * self.width_to_height_ratio
        Rectangle.__init__(self, **kwargs)
        if self.show_label:
            self.add_label()

    def add_label(self):
        label = TexMobject(self.label, color = self.label_color)
        label.scale(self.label_scale)
        label.move_to(self.get_center())
        self.add(label)


class Skabretti(Polygon):
    CONFIG = {
        "angle": 30,
        "show_angle_arc": True,
        "arc_radius": 1,
        "show_angle_label": True,
        "angle_label": "\\theta",
        "angle_label_scale": 1,
        "angle_label_buffer": MED_SMALL_BUFF,
        "kinetic_coefficient_of_friction": 0.5,
        "static_coefficient_of_friction": 0.5,
        "max_height": 5,
        "max_width": 11,
        "color": WHITE,
    }

    def __init__(self, **kwargs):
        digest_config(self, kwargs)

        tan = math.tan(self.angle * DEGREES)
        height, width = self.max_height, self.max_width
        if width * tan > height:
            width = height / tan
        else:
            height = width * tan

        vertices = [
            (0, 0, 0),
            (0, height, 0),
            (width, 0, 0)
        ]

        Polygon.__init__(self, *vertices)

        if self.show_angle_arc:
            self.add_angle_arc()
        if self.show_angle_label:
            self.add_angle_label()

    def add_angle_arc(self):
        ang = self.angle * DEGREES
        angle_arc = Arc(PI-ang, ang, arc_center = self.get_vertices()[2], radius = self.arc_radius)
        self.add(angle_arc)
    
    def add_angle_label(self):
        ang = self.angle * DEGREES
        angle_name = TexMobject(self.angle_label)
        angle_name.move_to(self.get_vertices()[2] + (1 + self.angle_label_buffer) * self.arc_radius * np.array((-math.cos(ang/2), math.sin(ang/2), 0)))
        #angle_name.scale(min(1, abs(angle_name.get_x())) * math.sin(theta))
        self.add(angle_name)


class KubburOgSkabretti(Scene):
    CONFIG = {
        # Breytur (í SI-einingum)
        "m": 1,
        "theta": 30,  # í gráðum
        "mu_k": 0.5,
        "mu_s": 0.5,
        "g": 9.82,

        # Stillingar
        "show_legend": True,
        "show_ind_force_vectors": False,
        "show_total_force_vector": False,
        "total_force_color": COLOR_MAP["RED_C"],
        "show_force_names": False,
    }

    def construct(self):
        skabretti = Skabretti()
        skabretti.to_corner(DL)

        ang = skabretti.angle*DEGREES
        kubbur = Kubbur()
        kubbur.rotate(-ang)
        kubbur.move_to(skabretti.get_vertices()[1])
        kubbur.shift(np.array([kubbur.height/2 * math.sin(ang), kubbur.height/2 * math.cos(ang), 0]))
        kubbur.shift(np.array([kubbur.width * math.cos(ang), -kubbur.width * math.sin(ang), 0]))

        if self.show_legend:
            self.add_legend()

        self.play(ShowCreation(skabretti), ShowCreation(kubbur))

    def add_legend(self, position=UR, scale=0.7,
                   show_m=True, show_theta=True, show_mu_k=True, show_mu_s=True, show_g=True):

        m_leg = TexMobject("m", "=", "\SI{%s}{kg}" % (self.m))
        m_leg[0].align_to(m_leg[1], RIGHT).next_to(m_leg[1], LEFT)
        m_leg[2].align_to(m_leg[1], LEFT).next_to(m_leg[1], RIGHT)

        theta_leg = TexMobject("\\theta", "=", "%s ^\\circ" % (self.theta))
        theta_leg[0].align_to(m_leg[0], RIGHT)
        theta_leg[1].align_to(m_leg[1], LEFT)
        theta_leg[2].align_to(m_leg[2], LEFT)

        mu_k_leg = TexMobject("\\mu_k", "=", "\\num{%s}" % (self.mu_k))
        mu_k_leg[0].align_to(m_leg[0], RIGHT)
        mu_k_leg[1].align_to(m_leg[1], LEFT)
        mu_k_leg[2].align_to(m_leg[2], LEFT)

        mu_s_leg = TexMobject("\\mu_s", "=", "\\num{%s}" % (self.mu_s))
        mu_s_leg[0].align_to(m_leg[0], RIGHT)
        mu_s_leg[1].align_to(m_leg[1], LEFT)
        mu_s_leg[2].align_to(m_leg[2], LEFT)

        g_leg = TexMobject("g", "=", "\\SI{%s}{\\unitfrac{m}{s^2}}" % (self.g))
        g_leg[0].align_to(m_leg[0], RIGHT)
        g_leg[1].align_to(m_leg[1], LEFT)
        g_leg[2].align_to(m_leg[2], LEFT)

        legend = VGroup()

        if show_m:
            legend.add(m_leg)
        if show_theta:
            legend.add(theta_leg)
        if show_mu_k:
            legend.add(mu_k_leg)
        if show_mu_s:
            legend.add(mu_s_leg)
        if show_g:
            legend.add(g_leg)

        self.add(legend.arrange(DOWN).scale(scale).to_corner(UR))


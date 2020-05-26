from manimlib.imports import *


class Block(Rectangle):
    CONFIG = {
        "mass": 1,
        "fill_opacity": 1.0,
        "fill_color": WHITE,
        "stroke_opacity": 0.0,
        "show_label": True,
        "label": "m",
        "max_label_scale": 1,
        "label_color": BLACK,
        "unit_height": 0.75,
        "width_to_height_ratio": (1 + 5 ** 0.5) / 2,
        "scale_with_mass": True,
    }

    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        self.height = self.unit_height
        if self.scale_with_mass:
            self.height *= self.mass ** 0.5
        self.width = self.height * self.width_to_height_ratio
        Rectangle.__init__(self, **kwargs)
        if self.show_label:
            self.add_label()

    def add_label(self):
        label = TexMobject(self.label, color=self.label_color)
        label.scale(min(self.max_label_scale, self.height))
        label.move_to(self.get_center())
        self.add(label)


class MovingBlock(Block):
    CONFIG = {
        "start_point": LEFT,
        "end_point": RIGHT,
        "align_block_edge_to_path": True,
        "distance_traveled": 0,
        "velocity": 0,
        "acceleration": 0,
    }

    def __init__(self, **kwargs):
        Block.__init__(self, **kwargs)
        path = self.get_def_path()
        ang = path.get_angle()
        l = path.get_length()

        if self.align_block_edge_to_path:
            path.shift(np.array([self.height/2 * math.sin(-ang), self.height/2 * math.cos(-ang), 0]))
            path.set_length(path.get_length()-self.width)
            l = path.get_length()

        self.rotate(ang)
        self.add_updater(lambda b: b.move_to(straight_path(
            path.get_start(), path.get_end(), self.distance_traveled / path.get_length())))
        
        self.true_path = path

    def get_def_path(self):
        return Line(self.start_point, self.end_point)

    def get_true_path(self):
        return self.true_path

    def move(self, d):
        self.distance_traveled = d


class InclinedPlane(Polygon):
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
        "max_height": 6,
        "max_width": 10,
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

        self.arc_radius = min(self.arc_radius, width/5)
        if self.show_angle_arc:
            self.add_angle_arc(width)
        if self.show_angle_label:
            self.add_angle_label()

    def add_angle_arc(self, width):
        ang = self.angle * DEGREES
        angle_arc = Arc(PI-ang, ang, arc_center=self.get_vertices()
                        [2], radius=self.arc_radius)
        self.add(angle_arc)

    def add_angle_label(self):
        ang = self.angle * DEGREES
        angle_name = TexMobject(self.angle_label)
        angle_name.move_to(self.get_vertices()[2] + (1 + self.angle_label_buffer)
                           * self.arc_radius * np.array((-math.cos(ang/2), math.sin(ang/2), 0)))
        angle_name.scale(min(self.angle_label_scale, abs(
            angle_name.get_x()-self.get_vertices()[2][0]) * math.tan(ang)))
        self.add(angle_name)


class BlockSlidingDownInclinedPlane(Scene):
    CONFIG = {
        # Values (in SI-units)
        "m": 1,
        "theta": 30,  # degrees
        "mu_k": 0.5,
        "mu_s": 0.5,
        "g": 9.82,

        "show_legend": True,
        "show_ind_force_vectors": False,
        "show_total_force_vector": False,
        "total_force_color": COLOR_MAP["RED_C"],
        "show_force_names": False,
    }

    def construct(self):
        incl_plane = InclinedPlane(
            angle=self.theta, kinetic_coefficient_of_friction=self.mu_k, static_coefficient_of_friction=self.mu_s)
        incl_plane.to_corner(DL)

        ang = incl_plane.angle*DEGREES

        block = MovingBlock(mass=self.m, start_point=incl_plane.get_vertices()[
                            1], end_point=incl_plane.get_vertices()[2])

        if self.show_legend:
            self.add_legend()

        self.add(incl_plane, block)
        self.wait(1)

        dist_tracker = ValueTracker(0)
        block.add_updater(lambda b: b.move(dist_tracker.get_value()))
        self.play(dist_tracker.set_value,
                  block.get_true_path().get_length(), run_time=3)
        self.wait(1)

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

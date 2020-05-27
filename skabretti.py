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
        "displacement": 0,
        "velocity": 0,
        "acceleration": 0,
    }

    def __init__(self, **kwargs):
        Block.__init__(self, **kwargs)
        path = Line(self.start_point, self.end_point)
        # self.displacement = self.initial_displacement
        ang = path.get_angle()
        self.path_length = path.get_length()

        if self.align_block_edge_to_path:
            path.shift(np.array([self.height/2 * math.sin(-ang), self.height/2 * math.cos(-ang), 0]))
            path.set_length(self.path_length-self.width)
            self.path_length = path.get_length()

        self.true_path = path
        self.rotate(ang)
        self.move_to(path.get_start())

        def update_block(self, dt):
            self.velocity += self.acceleration * dt
            self.displacement += self.velocity * dt
            self.move_to(straight_path(self.true_path.get_start(), self.true_path.get_end(), self.displacement / self.path_length))
            if self.displacement >= self.path_length:
                self.velocity = 0
                self.acceleration = 0
                self.move_to(self.true_path.get_end())

        self.add_updater(update_block)

    def start_move(self, a, v=0):
        self.acceleration = a
        self.velocity = v
    



class InclinedPlane(Polygon):
    CONFIG = {
        "angle": 30,
        "show_angle_arc": True,
        "arc_radius": 1,
        "show_angle_label": True,
        "angle_label": "\\theta",
        "angle_label_scale": 1,
        "angle_label_buffer": MED_SMALL_BUFF,
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
        "mu_k": 0.2,
        "mu_s": 0.2,
        "g": 9.82,
        "v0": 0,

        "combine_mu": True,
        "show_legend": True,
        "show_individual_force_vectors": True,
        "show_total_force_vector": True,
        "clip_forces_behind_block": False,
        "force_scale": 0.2,
        "hide_forces_in_slide": False,
        "show_force_names_by_vectors": True,
        "show_force_names_in_legend": False,
        "break_force_into_components": True,
        "force_colors": {
            "gravity": BLUE,
            "normal": YELLOW,
            "friction": RED,
            "total": GREEN,
        },
        "force_names": {
            "gravity": "\\vv{F_g}",
            "normal": "\\vv{Þ}",
            "friction": "\\vv{F_{\\text{nún}}}",
            "total": "\\vv{F_{\\text{heild}}}",
        },
        "time_before_slide": 1,
    }

    def construct(self):
        if self.combine_mu:
            self.mu_s = self.mu_k
        
        incl_plane = InclinedPlane(
            angle=self.theta, kinetic_coefficient_of_friction=self.mu_k, static_coefficient_of_friction=self.mu_s)
        incl_plane.to_corner(DL)

        ang = incl_plane.angle*DEGREES

        block = MovingBlock(mass=self.m, start_point=incl_plane.get_vertices()[
                            1], end_point=incl_plane.get_vertices()[2])
        
        a = self.g * (math.sin(ang) - self.mu_k * math.cos(ang))
        total_time = (-self.v0 + math.sqrt(self.v0**2 + 2*a*block.path_length)) / a
        if self.mu_s >= math.tan(ang):
            a = 0
            total_time = 0

        if self.show_legend:
            self.add_legend(self.combine_mu)

        self.add(incl_plane, block)
        if self.show_individual_force_vectors or self.show_total_force_vector:
            self.show_force_vectors(block)
        else:
            self.wait(self.time_before_slide)

        block.start_move(a, self.v0)
        self.wait(total_time + 1)

    def show_force_vectors(self, block):
        def get_polar_vector(length, angle):
            return Vector(length*np.array([math.cos(angle), math.sin(angle), 0]))

        grav = get_polar_vector(self.force_scale*self.m*self.g, -PI/2).set_color(self.force_colors['gravity'])
        norm = get_polar_vector(self.force_scale*self.m*self.g*math.cos(self.theta*DEGREES), PI/2-self.theta*DEGREES).set_color(self.force_colors['normal'])
        if self.mu_s >= math.tan(self.theta*DEGREES):
            fric = Vector(-(grav.get_end() + norm.get_end()))
        else:
            fric = get_polar_vector(self.force_scale*self.mu_k*self.m*self.g*math.cos(self.theta*DEGREES), PI-self.theta*DEGREES)
        fric.set_color(self.force_colors['friction'])
        
        forces = VGroup(grav, norm, fric)
        edges = VGroup()

        if self.show_total_force_vector:
            print("blinggg")
            total = get_polar_vector(self.force_scale*self.m*self.g*((math.sin(self.theta*DEGREES)-self.mu_k*math.cos(self.theta*DEGREES))), -self.theta*DEGREES)
            total.set_color(self.force_colors['total'])
            forces.add(total)

        if self.break_force_into_components:
            vert_comp = DashedLine(start=ORIGIN, end=-norm.get_end()).set_color(self.force_colors['gravity'])
            horz_comp = DashedLine(start=vert_comp.get_end(), end=grav.get_end()).set_color(self.force_colors['gravity'])
            forces.add(vert_comp, horz_comp)

            def update_horz(horz):
                horz.put_start_and_end_on(vert_comp.get_end(), grav.get_end())
            horz_comp.add_updater(update_horz)

        def move_to_block(vectors):
            for v in vectors:
                v.move_to(block.get_center()).shift(v.get_length()/2*v.get_unit_vector())

        forces.add_updater(move_to_block)

        if self.clip_forces_behind_block:
            self.add_foreground_mobject(block)
        self.play(ShowCreation(forces))
        self.wait(self.time_before_slide)
        if self.hide_forces_in_slide:
            self.remove(forces)

        

    def add_legend(self, combine_mu, position=UR, scale=0.7,
                   show_m=True, show_theta=True, show_mu_k=True, show_mu_s=True, show_g=False):

        m_leg = TexMobject("m", "=", "\SI{%s}{kg}" % (self.m))
        m_leg[0].align_to(m_leg[1], RIGHT).next_to(m_leg[1], LEFT)
        m_leg[2].align_to(m_leg[1], LEFT).next_to(m_leg[1], RIGHT)

        theta_leg = TexMobject("\\theta", "=", "%s ^\\circ" % (self.theta))
        theta_leg[0].align_to(m_leg[0], RIGHT)
        theta_leg[1].align_to(m_leg[1], LEFT)
        theta_leg[2].align_to(m_leg[2], LEFT)

        mu_leg = TexMobject("\\mu", "=", "\\num{%s}" % (self.mu_k))
        mu_leg[0].align_to(m_leg[0], RIGHT)
        mu_leg[1].align_to(m_leg[1], LEFT)
        mu_leg[2].align_to(m_leg[2], LEFT)
        
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
            if combine_mu:
                legend.add(mu_leg)
            else:
                legend.add(mu_k_leg)
        if show_mu_s and not combine_mu:
            legend.add(mu_s_leg)
        if show_g:
            legend.add(g_leg)

        self.add(legend.arrange(DOWN).scale(scale).to_corner(UR))

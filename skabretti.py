from manimlib.imports import *

### ATH ###
# Til þess að hraðagildið í legendinu líti ekki shitty út
# er hægt að breyta DecimalNumber klasanum (manimlib/mobject/numbers)
# á eftirfarandi hátt:
#
#     self.arrange(
#         buff=self.digit_to_digit_buff,
#         aligned_edge=DOWN
#     )
#
#     if self.unit is not None:
#         self.unit_sign = SingleStringTexMobject(self.unit, color=self.color)
#         self.unit_sign.next_to(self, RIGHT, buff=self.unit_buff) (*)
#         self.add(self.unit_sign)
#
# (sem sagt færa seinni blokkina undir fyrri og bæta við (*) í henni)
# Einnig þarf að bæta '"unit_buff": 0.15,' við CONFIG dictið.

class Block(Rectangle):
    CONFIG = {
        "mass": 1,
        "fill_opacity": 1.0,
        "fill_color": WHITE,
        "stroke_opacity": 0.0,
        "show_label": True,
        "label_text": "m",
        "label_scale": 1,
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
        self.label = self.get_label()

    def get_label(self):
        return TexMobject(self.label_text, color=self.label_color) \
               .scale(self.label_scale * self.height) \
               .add_updater(lambda l: l.move_to(self.get_center()))

    def add_label(self, scene, rotate_by=0, animate=False, animation=Write):
        self.label.rotate(rotate_by)
        if animate:
            scene.play(animation(self.label))
        else:
            scene.add(self.label)

    def remove_label(self, scene, animate=False, animation=Uncreate):
        if animate:
            scene.play(animation(self.label))
        else:
            scene.remove(self.label)

    # ew
    def get_critical_point_on_block(self, direction):
        if (direction == ORIGIN).all():
            return (self.get_vertices()[0] + self.get_vertices()[1] + self.get_vertices()[2] + self.get_vertices()[3]) / 4
        if (direction == RIGHT).all():
            return mid(self.get_vertices()[1], self.get_vertices()[2])
        if (direction == UR).all():
            return self.get_vertices()[1]
        if (direction == UP).all():
            return mid(self.get_vertices()[0], self.get_vertices()[1])
        if (direction == UL).all():
            return self.get_vertices()[0]
        if (direction == LEFT).all():
            return mid(self.get_vertices()[0], self.get_vertices()[3])
        if (direction == DL).all():
            return self.get_vertices()[3]
        if (direction == DOWN).all():
            return mid(self.get_vertices()[2], self.get_vertices()[3])
        if (direction == DR).all():
            return self.get_vertices()[2]
        else:
            raise Exception("Invalid input")


class MovingBlock(Block):
    CONFIG = {
        "start_point": LEFT,
        "end_point": RIGHT,
        "align_block_edge_to_path": True,
        "displacement": 0,
        "velocity": 0,
        "acceleration": 0,
    }

    FORCES = {
        "LABELS": {},
        "COMPS": {}
    }

    def __init__(self, **kwargs):
        Block.__init__(self, **kwargs)
        path = Line(self.start_point, self.end_point)
        ang = path.get_angle()
        self.path_length = path.get_length()

        if self.align_block_edge_to_path:
            path.shift(
                np.array([self.height/2 * math.sin(-ang), self.height/2 * math.cos(-ang), 0]))
            path.set_length(self.path_length-self.width)
            self.path_length = path.get_length()

        self.true_path = path
        self.rotate(ang)
        self.move_to(path.get_start())

        self.stopped, self.last_velocity = False, 0
        def update_block(self, dt):
            if self.displacement >= self.path_length:
                self.last_velocity = self.velocity
                self.velocity = 0
                self.acceleration = 0
                self.displacement = self.path_length

            self.velocity += self.acceleration * dt
            self.displacement += self.velocity * dt

            self.move_to(straight_path(
                    self.true_path.get_start(), 
                    self.true_path.get_end(), 
                    self.displacement / self.path_length
                ))

        self.add_updater(update_block)

    def start_move(self, a, v=0):
        self.acceleration = a
        self.velocity = v

    def move_to_alpha(self, alpha):
        self.displacement = alpha * self.path_length

    def make_force_vector(self, name, vector, position=ORIGIN, **kwargs):
        vec = Vector(**kwargs)
        vec.updater = lambda v: v.put_start_and_end_on(
                self.get_critical_point_on_block(position), 
                self.get_critical_point_on_block(position) + vector)
        vec.updater(vec)
        self.FORCES[name] = vec

    def add_force_vectors(self, scene, *names, animate=True, animation=ShowCreation, update=True):
        vecs = [self.FORCES[name] for name in names]
        for vec in vecs:
            if update:
                vec.add_updater(vec.updater)
        if animate:
            scene.play(*[animation(vec) for vec in vecs])
        else:
            scene.add(*vecs)
        # for vec in vecs:
        #     vec.resume_updating()

    def remove_force_vectors(self, scene, *names, animate=True, animation=Uncreate):
        vecs = [self.FORCES[name] for name in names]
        for vec in vecs:
            if vec.updater in vec.get_updaters():
                vec.remove_updater(vec.updater)
        if animate:
            scene.play(*[animation(vec) for vec in vecs])
        else:
            scene.remove(*vecs)

    def make_force_label(self, name, text, direction, scale, **kwargs):
        label = TexMobject(text, **kwargs).scale(scale)
        label.updater = lambda l: l.next_to(self.FORCES[name].get_end(), direction)
        label.updater(label)
        self.FORCES["LABELS"][name] = label

    def add_force_labels(self, scene, *names, animate=True, animation=Write, update=True):
        labels = [self.FORCES["LABELS"][name] for name in names]
        for label in labels:
            if update:
                label.add_updater(label.updater)
        if animate:
            scene.play(*[animation(label) for label in labels])
        else:
            scene.add(*labels)

    def remove_force_labels(self, scene, *names, animate=True, animation=Uncreate):
        labels = [self.FORCES["LABELS"][name] for name in names] 
        for label in labels:
            if label.updater in label.get_updaters():
                label.remove_updater(label.updater)       
        if animate:
            scene.play(*[animation(label) for label in labels])
        else:
            scene.remove(*labels)
        

class InclinedPlane(Polygon):
    CONFIG = {
        "angle": 30,
        "show_angle_arc": True,
        "arc_radius": 1,
        "show_angle_label": True,
        "angle_label_text": "\\theta",
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
        # if self.show_angle_arc:
        #     self.add(self.angle_arc)
        # if self.show_angle_label:
        #     self.add(self.angle_label)

    def get_angle_arc(self):
        return Arc(PI-self.angle*DEGREES, self.angle*DEGREES, arc_center=self.get_vertices()[2], radius=self.arc_radius)

    def get_angle_label(self):
        ang = self.angle * DEGREES
        angle_label = TexMobject(self.angle_label_text)
        angle_label.move_to(self.get_vertices()[2] + (1 + self.angle_label_buffer) * self.arc_radius * np.array((-math.cos(ang/2), math.sin(ang/2), 0)))
        angle_label.scale(min(self.angle_label_scale, abs(angle_label.get_x()-self.get_vertices()[2][0]) * math.tan(ang)))
        return angle_label

    def add_angle_arc_and_label(self, scene, animate=True, arc_animation=ShowCreation, label_animation=Write):
        self.angle_arc = self.get_angle_arc()
        self.angle_label = self.get_angle_label()
        if animate:
            scene.play(AnimationGroup(
                arc_animation(self.angle_arc), 
                label_animation(self.angle_label), 
                lag_ratio=0.6
            ))
        else:
            scene.add(self.angle_arc, self.angle_label)
    
    def remove_angle_arc_and_label(self, scene, animate=True, animation=Uncreate):
        if animate:
            scene.play(AnimationGroup(
                animation(self.angle_label), 
                animation(self.angle_arc), 
                lag_ratio=0.6
            ))
        else:
            scene.remove(self.angle_arc, self.angle_label)

    def set_angle(self, angle):
        new_height = self.get_width()*math.tan(angle)
        self.get_vertices()[1]



class BlockAndInclinedPlaneScene(Scene):
    CONFIG = {
        # Values (in SI-units)
        "m": 1,
        "theta": 30, # degrees
        "mu_k": 0.2,
        "mu_s": 0.5,
        "g": 9.82,
        "v0": 0,
        "s0": 0,

        "combine_mu": True,
        "show_legend": True,
        "show_in_legend": {
            "m": False,
            "theta": False,
            "mu_k": True,
            "mu_s": False,
            "g": False,
            "av_seperator": True,
            "a": True,
            "v": True,
            "force_seperator": False,
            "forces": False,
        },
        "legend_scale": 0.8,
        # "show_individual_force_vectors": False,
        # "show_total_force_vector": False,
        "show_forces": {
            "gravity": True,
            "normal": True,
            "friction": True,
            "total": True,
        },
        "set_forces_to_block_edge": True,
        "force_scale": 0.22,
        "hide_forces_in_slide": True,
        "show_force_names_by_vectors": True,
        "break_force_into_components": True,
        "force_colors": {
            "gravity": BLUE,
            "normal": YELLOW,
            "friction": RED,
            "total": GREEN,
        },
        "force_names": {
            "gravity": "F_g",
            "normal": "\\text{\\textit{Þ}}",
            "friction": "F_{\\text{nún}}",
            "total": "F_{\\text{heild}}",
        },
        "time_before_slide": 1,
    }

    def construct(self):
        if self.combine_mu:
            self.mu_s = self.mu_k

        incl_plane = InclinedPlane(angle=self.theta)
        incl_plane.to_corner(DL)

        ang = self.theta*DEGREES

        block = MovingBlock(
            mass=self.m, 
            start_point=incl_plane.get_vertices()[1], 
            end_point=incl_plane.get_vertices()[2], 
            displacement=self.s0
        )

        self.stopped = self.mu_s >= math.tan(ang)

        if self.stopped:
            a = 0
            total_time = 0
        else:
            a = self.g * (math.sin(ang) - self.mu_k * math.cos(ang))
            total_time = (-self.v0 + math.sqrt(self.v0 ** 2 + 2*a*block.path_length)) / a

        if self.show_legend:
            self.add(self.get_legend(block))

        self.add(incl_plane, block)
        # if self.show_individual_force_vectors or self.show_total_force_vector:  
        self.show_force_vectors_new(block)
        # else:
        #     self.wait(self.time_before_slide)

        block.start_move(a, self.v0)
        self.wait(total_time + 1)

    # Frekar ljótur kóði en whatever
    # def show_force_vectors(self, block):  
    #     ang = self.theta * DEGREES
    #     mg = self.force_scale * self.m * self.g
    #     mgsin = mg * math.sin(ang)
    #     mgcos = mg * math.cos(ang)

    #     def get_polar_vector(length, angle):
    #         return Vector(length*np.array([math.cos(angle), math.sin(angle), 0]))

    #     grav = get_polar_vector(mg, -PI/2).set_color(self.force_colors['gravity'])
    #     norm = get_polar_vector(mgcos, PI/2-ang).set_color(self.force_colors['normal'])
    #     if self.mu_s >= math.tan(ang):
    #         fric = get_polar_vector(mgsin, PI-ang)
    #         total = Vector(ORIGIN)
    #     else:
    #         fric = get_polar_vector(self.mu_k*mgcos, PI-ang)
    #         total = get_polar_vector(
    #             mg*((math.sin(ang)-self.mu_k*math.cos(ang))), -ang)
    #     fric.set_color(self.force_colors['friction'])
    #     total.set_color(self.force_colors['total'])

    #     forces = VGroup()
    #     if self.show_individual_force_vectors:
    #         forces.add(grav, norm, fric)
    #     if self.show_total_force_vector:
    #         forces.add(total)

    #     vert_comp = None
    #     horz_comp = None
    #     if self.break_force_into_components and self.show_individual_force_vectors:
    #         horz_comp = DashedLine(start=vert_comp.get_end(), end=grav.get_end(), stroke_width=2).set_color(self.force_colors['gravity'])
    #         forces.add(vert_comp, horz_comp)

    #         def update_horz(horz):
    #             horz.put_start_and_end_on(vert_comp.get_end(), grav.get_end())
    #         horz_comp.add_updater(update_horz)

    #     # ég hef ekki hugmynd af hverju það er 2.33 og 3.33 hérna, ætti bæði að vera 2 ¯\_(ツ)_/¯
    #     bottom_offset = -block.get_height()/3.33 * \
    #         np.array([math.sin(ang), math.cos(ang), 0])
    #     top_offset = -bottom_offset
    #     left_offset = -block.get_width()/2.33 * \
    #         np.array([math.cos(-ang), math.sin(-ang), 0])
    #     right_offset = -left_offset

    #     grav_name = TexMobject("\\vec{", self.force_names['gravity']).set_color(self.force_colors['gravity'], "}").scale(min(1, grav.get_length()))
    #     norm_name = TexMobject("\\vec{", self.force_names['normal']).set_color(self.force_colors['normal'], "}").scale(min(1, norm.get_length()))
    #     fric_name = TexMobject("\\vec{", self.force_names['friction']).set_color(self.force_colors['friction'], "}").scale(min(1, fric.get_length()))
    #     total_name = TexMobject("\\vec{", self.force_names['total']).set_color(self.force_colors['total'], "}").scale(min(1, total.get_length()))
    #     names = VGroup()
    #     if self.show_individual_force_vectors:
    #         names.add(grav_name, norm_name, fric_name)
    #     if self.show_total_force_vector:
    #         names.add(total_name)
    #     force_names = {
    #         grav: grav_name,
    #         norm: norm_name,
    #         fric: fric_name,
    #         total: total_name,
    #     }

    #     def move_to_block(vectors):
    #         force_positions = {
    #             grav: block.get_center(),
    #             norm: block.get_center(),
    #             fric: block.get_center(),
    #             total: block.get_center(),
    #             vert_comp: block.get_center(),
    #             horz_comp: block.get_center(),
    #         }
    #         if self.set_forces_to_block_edge:
    #             force_positions[grav] += bottom_offset
    #             force_positions[norm] += top_offset
    #             force_positions[fric] += left_offset
    #             force_positions[total] += right_offset
    #             force_positions[vert_comp] += bottom_offset

    #         for v in vectors:
    #             v.move_to(force_positions[v]).shift(
    #                 v.get_length()/2*v.get_unit_vector())
    #             if self.show_force_names_by_vectors and v not in (vert_comp, horz_comp):
    #                 side = RIGHT
    #                 if v in (fric):
    #                     side = UP
    #                 force_names[v].next_to(v.get_end(), side)

    #     move_to_block(forces)
    #     if self.set_forces_to_block_edge:
    #         self.add_foreground_mobject(block)

    #     self.play(ShowCreation(forces))
    #     if self.show_force_names_by_vectors:
    #         self.play(Write(names))
    #     forces.add_updater(move_to_block)
    #     self.wait(self.time_before_slide)

    #     if self.hide_forces_in_slide:
    #         forces.remove_updater(move_to_block)
    #         if self.show_force_names_by_vectors:
    #             self.play(Uncreate(names))
    #         self.play(Uncreate(forces))

    def show_force_vectors_new(self, block):
        ang = self.theta * DEGREES
        mg = self.force_scale * self.m * self.g
        mgsin = mg * math.sin(ang)
        mgcos = mg * math.cos(ang)

        force_sizes = {
            "gravity": mg,
            "normal": mgcos,
            "friction": mgsin if self.stopped else self.mu_k*mgcos,
            "total": 0 if self.stopped else mg*((math.sin(ang)-self.mu_k*math.cos(ang)))
        }

        force_directions = {
            "gravity": DOWN,
            "normal": np.array([math.cos(PI/2 - ang), math.sin(PI/2 - ang), 0]),
            "friction": np.array([math.cos(PI - ang), math.sin(PI - ang), 0]),
            "total": np.array([math.cos(-ang), math.sin(-ang), 0])
        }

        force_origins = {
            "gravity": DOWN,
            "normal": UP,
            "friction": LEFT,
            "total": RIGHT
        }

        for f in ("gravity", "normal", "friction", "total"):
            block.make_force_vector(f, 
                                    force_sizes[f] * force_directions[f], 
                                    force_origins[f] if self.set_forces_to_block_edge else ORIGIN, 
                                    color=self.force_colors[f])

        block.add_force_vectors(self, *[k for k,v in self.show_forces.items() if v], animate=True)

        self.wait(self.time_before_slide)

        if self.hide_forces_in_slide:
            block.remove_force_vectors(self, *[k for k,v in self.show_forces.items() if v], animate=False)
            self.wait(1)
            block.add_force_vectors(self, *[k for k,v in self.show_forces.items() if v], animate=True)
            


    def get_legend(self, block):
        mg = self.m * self.g
        mgsin = mg * math.sin(self.theta*DEGREES)
        mgcos = mg * math.cos(self.theta*DEGREES)
        a = round(self.g * (math.sin(self.theta*DEGREES) - self.mu_k * math.cos(self.theta*DEGREES)), 2)

        vals = {
            "m": TexMobject("\SI{%s}{kg}" % (self.m)),
            "theta": TexMobject("\\num{%s}^\\circ" % (self.theta)),
            "mu_k": TexMobject("\\num{%s}" % (self.mu_k)),
            "mu_s": TexMobject("\\num{%s}" % (self.mu_s)),
            "g": TexMobject("\\SI{%s}{\\unitfrac{m}{s^2}}" % (self.g)),
            "a": TexMobject("\\SI{%s}{\\unitfrac{m}{s^2}}" % (max(0, a))),
            "v": DecimalNumber(block.velocity, unit="\\unitfrac{m}{s}}", digit_to_digit_buff=0.06)
                 .add_updater(lambda v: v.set_value(block.velocity if not block.stopped else block.last_velocity)),
            "forces": {
                "gravity": TexMobject("\\SI{%s}{N}" % (round(mg, 2))),
                "normal": TexMobject("\\SI{%s}{N}" % (round(mgcos, 2))),
                "friction": TexMobject("\\SI{%s}{N}" % (round(mgsin, 2))),
                "total": TexMobject("\\SI{%s}{N}" % (max(0, round(self.m*a, 2))))
            }
        }
        if self.stopped:
            vals["v"] = TexMobject("\\SI{0}{\\unitfrac{m}{s}}")

        last = ORIGIN
        legend = VGroup()

        def add_to_legend(c_leg, mylast):
            c_leg[1].next_to(mylast, DOWN, buff=MED_LARGE_BUFF)
            c_leg[0].align_to(c_leg[1], RIGHT).next_to(c_leg[1], LEFT)
            c_leg[2].align_to(c_leg[1], LEFT).next_to(c_leg[1], RIGHT)
            legend.add(c_leg)
            mylast = c_leg[1].get_center()
            return mylast

        for comp in self.show_in_legend:
            if self.show_in_legend[comp]:
                if comp == "forces":
                    for f_comp in vals[comp]:
                        last = add_to_legend(TexMobject(self.force_names[f_comp], "=").add(vals[comp][f_comp]).set_color(self.force_colors[f_comp]), last)
                    continue

                if comp in ("av_seperator", "force_seperator"):
                    sep = Line(legend.get_left(), legend.get_right(), stroke_width=1)
                    sep.next_to(last, DOWN, buff=MED_LARGE_BUFF).shift(np.array([legend.get_x()-sep.get_x(), 0, 0]))
                    sep.add_updater(lambda s: s.set_width(legend.get_width()).shift(np.array([legend.get_center()[0]-s.get_center()[0], 0, 0])))
                    legend.add(sep)
                    last = last + np.array([0, sep.get_y()-last[1], 0])
                    continue

                name = comp
                if comp == "mu_k" and self.combine_mu:
                    name = "\\mu"
                elif comp in ("theta", "mu_k", "mu_s"):
                    name = "\\" + name
                    self.show_in_legend["mu_s"] = False 
                
                last = add_to_legend(TexMobject(name, "=").add(vals[comp]), last)

        return legend.scale(self.legend_scale).to_corner(UR)


class BlockInclinedPlaneDerivation(MovingCameraScene):
    CONFIG = {
        "m": 1,
        "theta": 30,
        "mu": 0.25,
        "g": 9.82,
        "force_scale": 0.2
    }

    def construct(self):
        plane = InclinedPlane(
            angle=self.theta, 
            show_angle_arc=False, 
            show_angle_label=False,
            arc_radius=1.5
        ).to_corner(DL)

        block = MovingBlock(
            mass=self.m, 
            start_point=plane.get_vertices()[1], 
            end_point=plane.get_vertices()[2], 
            displacement=1,
            label_scale=1.25,
            unit_height=0.9
        )

        ang = self.theta*DEGREES
        mg = self.force_scale * self.m * self.g
        mgcos = mg * math.cos(ang)

        self.forces = {
            "gravity": {
                "size": mg,
                "angle": -PI/2,
                "origin": DOWN,
                "color": BLUE,
                "label_text": "F_g",
                "label_direction": RIGHT
            },
            "normal": {
                "size": mgcos,
                "angle": PI/2 - ang,
                "origin": UP,
                "color": YELLOW,
                "label_text": "\\text{\\textit{Þ}}",
                "label_direction": RIGHT
            },
            "friction": {
                "size": self.mu*mgcos,
                "angle": PI - ang,
                "origin": LEFT,
                "color": RED,
                "label_text": "F_{\\text{nún}}",
                "label_direction": UL
            },
            "total": {
                "size": mg*((math.sin(ang)-self.mu*math.cos(ang))),
                "angle": -ang,
                "origin": RIGHT,
                "color": GREEN,
                "label_text": "F_{\\text{heild}}",
                "label_direction": RIGHT
            }
        }

        self.add(plane, block)

        self.add_sound("skabretti1.m4a")
        self.wait(4.75)

        block.add_label(self, rotate_by=-ang, animate=True)
        self.add_foreground_mobjects(block, block.label)
        self.wait(1.75)
        plane.add_angle_arc_and_label(self)
        self.wait(1.25)

        a = self.g * (math.sin(ang) - self.mu * math.cos(ang))
        slide_time = math.sqrt(2*block.path_length/a)

        block.start_move(a)
        self.wait(slide_time+2.5)

        self.add_sound("skabretti2.m4a")
        self.play(ApplyMethod(block.move_to_alpha, 0.5))
        block.move_to_alpha(0.5) # veit ekki afh þetta þarf
        self.wait(2.5)

        plane.remove_angle_arc_and_label(self)
        self.wait(0.75)

        self.play(
            self.camera_frame.set_width, plane.get_width(),
            self.camera_frame.move_to, block.get_center(),
            run_time=1.5
        )
        self.wait(2)

        self.make_forces(block)

        for f, sound in zip(("gravity", "normal", "friction"), ("skabretti3.m4a", "skabretti4.m4a", "skabretti5.m4a")):
            self.add_sound(sound)
            self.wait(2 if f=="gravity" else 1.25)
            block.add_force_vectors(self, f, update=False)
            block.add_force_labels(self, f, update=False)
            if f in ("gravity", "friction"):
                eq_text = {
                    "gravity": "mg",
                    "friction": "\\mu\\text{\\textit{Þ}}"
                }
                eq_mob = TexMobject("=", eq_text[f], color=self.forces[f]["color"]) \
                     .scale(min(0.8, self.forces[f]["size"]), about_point=block.FORCES["LABELS"][f].get_right()) \
                     .next_to(block.FORCES["LABELS"][f], RIGHT, aligned_edge=DOWN, buff=min(0.18, 0.18*self.forces[f]["size"]))
                self.wait(0.5 if f=="gravity" else 2.5)
                self.play(Write(eq_mob))
                self.wait(4 if f=="gravity" else 10)
                self.play(Uncreate(eq_mob))
            else:
                self.wait(3)

        self.wait(3)
        self.add_sound("skabretti6.m4a")
        self.wait(0.5)
        self.play(self.camera_frame.shift, 2.5*RIGHT, run_time=1.25)
        
        total_force_eq = TexMobject(
            self.forces["total"]["label_text"], 
            "=", self.forces["gravity"]["label_text"], 
            "+", self.forces["normal"]["label_text"], 
            "+", self.forces["friction"]["label_text"]).scale(0.75)
        total_force_eq.move_to(self.camera_frame.get_center() + np.array([
            self.camera_frame.get_width()/2 - 1,
            self.camera_frame.get_height()/2 - 1,
            0
        ]), aligned_edge=UR)
        for i, f in zip((0, 2, 4, 6), ("total", "gravity", "normal", "friction")):
            total_force_eq[i].set_color(self.forces[f]["color"])
        self.play(AnimationGroup(
            Write(total_force_eq[0:2]),
            AnimationGroup(
                ReplacementTransform(block.FORCES["LABELS"]["gravity"], total_force_eq[2], run_time=1.5),
                ReplacementTransform(block.FORCES["LABELS"]["normal"], total_force_eq[4], run_time=1.5),
                ReplacementTransform(block.FORCES["LABELS"]["friction"], total_force_eq[6], run_time=1.5),
                *[FadeIn(total_force_eq[i]) for i in (3, 5)],
                lag_ratio=0.2
            ),
            lag_ratio=0.2
        ))
        self.wait(1)

        vec_arrows = {}
        for i, f in zip((0, 2, 4, 6), ("total", "gravity", "normal", "friction")):
            arrow = Vector(tip_length=0.045).scale(0.25).set_color(self.forces[f]["color"]).next_to(total_force_eq[i], UP, aligned_edge=LEFT, buff=0.075).shift(0.04*RIGHT)
            vec_arrows[f] = arrow
        self.play(AnimationGroup(*[ShowCreation(arr) for arr in vec_arrows.values()], lag_ratio=0.2))
        self.wait(0.2)

        for f in ("gravity", "normal", "friction"):
            block.FORCES[f].generate_target()
        block.FORCES["gravity"].target.move_to(block.get_critical_point_on_block(RIGHT)).shift(self.forces["gravity"]["size"]/2*DOWN)
        block.FORCES["normal"].target.move_to(block.FORCES["gravity"].target.get_end(), aligned_edge=DL)
        block.FORCES["friction"].target.move_to(block.FORCES["normal"].target.get_end(), aligned_edge=DR)
        self.play(AnimationGroup(
            *[MoveToTarget(block.FORCES[f]) for f in ("gravity", "normal", "friction")],
            lag_ratio=0.4
        ))
        self.wait(1)

        block.add_force_vectors(self, "total", update=False)
        self.play(AnimationGroup(
            *[FadeOut(block.FORCES[f]) for f in ("gravity", "normal", "friction")],
            ReplacementTransform(total_force_eq[0], block.FORCES["LABELS"]["total"]),
            *[FadeOut(total_force_eq[i]) for i in range(1, 7)],
            *[FadeOut(arr) for arr in vec_arrows.values()]
        ))
        self.wait(2.5)

        self.add_sound("skabretti7.m4a")
        self.wait(0.75)

        self.camera_frame.save_state()
        self.play(
            self.camera_frame.move_to, block.FORCES["total"],
            self.camera_frame.set_width, block.FORCES["total"].get_width()*13
        )
        self.wait(1)

        total_force_projection_corner = np.array([block.FORCES["total"].get_end()[0], block.FORCES["total"].get_start()[1], 0])
        total_force_vert_comp = DashedLine(total_force_projection_corner, block.FORCES["total"].get_end()).set_color(self.forces["total"]["color"])
        total_force_horz_comp = DashedLine(block.FORCES["total"].get_start(), total_force_projection_corner).set_color(self.forces["total"]["color"])
        self.play(ShowCreation(total_force_horz_comp))
        self.wait(0.15)
        self.play(ShowCreation(total_force_vert_comp))
        self.wait(1)

        self.play(ApplyMethod(self.camera_frame.restore))

        total_force_vert_comp.add_updater(lambda c: c.put_start_and_end_on(
            np.array([block.FORCES["total"].get_end()[0], block.FORCES["total"].get_start()[1], 0]),
            block.FORCES["total"].get_end()
        ))
        total_force_horz_comp.add_updater(lambda c: c.put_start_and_end_on(
            block.FORCES["total"].get_start(),
            np.array([block.FORCES["total"].get_end()[0], block.FORCES["total"].get_start()[1], 0])
        ))

        x_vec = Vector(RIGHT, tip_length=0.15, stroke_width=1.25)
        y_vec = x_vec.copy().rotate(PI/2, about_point=ORIGIN)
        x_label = TexMobject("x").scale(0.6).next_to(x_vec, RIGHT, buff=SMALL_BUFF)
        y_label = TexMobject("y").scale(0.6).next_to(y_vec, UP, buff=SMALL_BUFF)
        block.remove_force_labels(self, "total", animation=FadeOut)
        axes = VGroup(x_vec, y_vec, x_label, y_label).scale(0.7).move_to(block.get_center()+3*RIGHT)
        self.play(ShowCreation(axes))
        self.wait(2)
        self.play(Rotate(axes, -ang, about_point=x_vec.get_start()), run_time=1.5)
        self.wait(2)

        plane.scale(4, about_point=block.get_critical_point_on_block(DOWN))

        self.play(
            *[Rotate(mob, ang, about_point=block.get_center()) for mob in 
                (plane, block, block.label, block.FORCES["total"], axes)],
            run_time=2.5
        )
        total_force_vert_comp.clear_updaters()
        total_force_horz_comp.clear_updaters()

        self.add_sound("skabretti8.m4a")
        self.wait(1)
        self.play(*[FadeOut(thing) for thing in (total_force_horz_comp, total_force_vert_comp)])
        self.wait(4)

        self.remove(block, block.label)
        block.remove_force_vectors(self, "total", animate=False)
        block = MovingBlock(
            mass=self.m, 
            start_point=plane.get_vertices()[1], 
            end_point=plane.get_vertices()[2], 
            label_scale=1.25,
            unit_height=0.9
        )
        block.move_to_alpha(0.5)
        block.add_label(self)
        for f in self.forces.keys():
            self.forces[f]["angle"] += ang
        self.make_forces(block)
        block.add_force_vectors(self, "total", animate=False, update=True)
        self.add(block, block.label)

        block.start_move(a)
        self.wait((math.sqrt(2*a*(block.path_length-block.displacement))) / a + 1.5)
        self.remove(block, block.FORCES["total"], block.label)
        block.move_to_alpha(0.5)
        self.add_sound("skabretti9.m4a")
        self.play(
            FadeIn(block),
            FadeIn(block.FORCES["total"]),
            FadeIn(block.label)
        )
        block.move_to_alpha(0.5)
        self.wait(4)

        self.add_sound("2logmal.m4a")
        self.wait(2)

        self.play(
            FadeOut(axes),
            self.camera_frame.set_width, self.camera_frame.get_width() * 1.3,
            self.camera_frame.shift, 1.5*RIGHT+UP
        )
        self.wait(2)

        newt_eq = TexMobject("F_{\\text{heild}}", "=", "m", "a")
        newt_eq.move_to(self.camera_frame.get_center() + np.array([
            self.camera_frame.get_width()/2 - 1.5,
            self.camera_frame.get_height()/2 - 1.5,
            0
        ]), aligned_edge=UR)
        newt_eq[0].set_color(self.forces["total"]["color"])
        vec_arrows = {}
        for i, f in zip((0, 3), ("total", "a")):
            arrow = Vector(tip_length=0.07).scale(0.35-0.05*i).set_color(self.forces[f]["color"] if f in self.forces.keys() else WHITE) \
                    .next_to(newt_eq[i], UP, aligned_edge=LEFT, buff=0.075).shift((0.04-0.01*i)*RIGHT)
            vec_arrows[f] = arrow
        self.play(AnimationGroup(
            ReplacementTransform(block.FORCES["total"].copy().clear_updaters(), vec_arrows["total"]),
            FadeIn(newt_eq[0]),
            ShowCreation(vec_arrows["a"]),
            Write(newt_eq[1:]),
            lag_ratio=0.7
        ))
        self.wait(7.5)

        new_force_eq = total_force_eq[2:].scale(1/0.75).next_to(newt_eq[1], LEFT)
        for i, f in zip((0, 2, 4), ("gravity", "normal", "friction")):
            arrow = Vector(tip_length=0.07).scale(0.35).set_color(self.forces[f]["color"]) \
                    .next_to(new_force_eq[i], UP, aligned_edge=LEFT, buff=0.075).shift((0.04)*RIGHT)
            vec_arrows[f] = arrow
        total_force = block.FORCES["total"]
        self.make_forces(block)
        self.play(AnimationGroup(
            ReplacementTransform(VGroup(newt_eq[0], vec_arrows["total"]), new_force_eq),
            Uncreate(total_force.clear_updaters()),
            *[ShowCreation(block.FORCES[f]) for f in ("gravity", "normal", "friction")],
            *[ShowCreation(arrow) for arrow in [v for k,v in vec_arrows.items() if k not in ("a", "total")]],
            lag_ratio=0.05
        ))
        self.wait(5)

        projection = block.FORCES["gravity"].copy()
        projection.rotate(-ang, about_point=projection.get_start()).scale(math.cos(-ang), about_point=projection.get_start())
        vert_comp = DashedLine(start=projection.get_start(), end=projection.get_end())
        horz_comp = DashedLine(start=projection.get_end(), end=block.FORCES["gravity"].get_end())
        comps = VGroup(vert_comp, horz_comp).set_color(self.forces["gravity"]["color"])
        self.play(ShowCreation(comps))
        self.wait(1)

        self.add_sound("skabretti10.m4a")
        self.wait(4)

        self.camera_frame.save_state()
        removed_in_grav_section = [r for r in self.get_mobjects() if r not in (block.FORCES["gravity"], comps, plane, self.camera_frame)]
        self.play(
            *[FadeOut(thing) for thing in removed_in_grav_section],
            self.camera_frame.move_to, block.FORCES["gravity"].get_center(),
            self.camera_frame.set_width, self.camera_frame.get_width()*0.7,
            lag_ratio=0.15,
            run_time=2
        )
        self.wait(1)
        ang_line = Line(LEFT, RIGHT)
        ang_line.set_angle(ang+PI)
        ang_line.set_length(self.forces["gravity"]["size"]/math.tan(ang)) 
        ang_line.move_to(block.FORCES["gravity"].get_end(), aligned_edge=DL)
        ang_line_arc = Arc(PI, ang, arc_center=ang_line.get_start(), radius=0.6)
        ang_line_label = TexMobject("\\theta").scale(0.6).move_to(ang_line.get_start()).shift(0.85*np.array([math.cos(PI+ang/2), math.sin(PI+ang/2), 0]))
        grav_line = Line(block.FORCES["gravity"].get_start(), block.FORCES["gravity"].get_end()).match_style(block.FORCES["gravity"])
        self.add_foreground_mobject(grav_line)
        self.play(AnimationGroup(
            ApplyMethod(self.camera_frame.move_to, VGroup(comps, block.FORCES["gravity"], ang_line, ang_line_arc).get_center()),
            ShowCreation(ang_line), 
            ShowCreation(ang_line_arc),
            FadeIn(ang_line_label),
            FadeOut(block.FORCES["gravity"].get_tip()),
            lag_ratio=0.15,
            run_time=2
        ))
        self.wait(1)

        right_angle_plane = Elbow(width=1).scale(0.2).move_to(ang_line.get_end(), aligned_edge=DL).rotate(ang, about_point=ang_line.get_end())
        right_angle_grav = Elbow(width=1, color=self.forces["gravity"]["color"]).scale(0.15).move_to(vert_comp.get_end(), aligned_edge=DL)
        similar_angle_plane = VGroup(*[Arc(0, ang-PI/2, arc_center=grav_line.get_start(), radius=r, stroke_width=2) for r in (0.3, 0.25)])
        similar_angle_grav = VGroup(*[Arc(PI, ang-PI/2, arc_center=grav_line.get_end(), radius=r, stroke_width=2, color=self.forces["gravity"]["color"]) for r in (0.25, 0.2)])

        self.play(AnimationGroup(
            ShowCreation(right_angle_plane),
            ShowCreation(right_angle_grav),
            ShowCreation(similar_angle_plane),
            ShowCreation(similar_angle_grav),
            lag_ratio=2,
            run_time=5
        ))
        self.wait(2)
        
        ang_line_arc_copy = ang_line_arc.copy()
        moving_triangle = VGroup(
            ang_line.copy(), 
            ang_line_arc_copy, 
            grav_line.copy(), 
            Line(ang_line.get_start(), grav_line.get_start()),
        )
        moving_triangle.stroke_opacity = 0.2

        moving_triangle.generate_target()
        moving_triangle.target.rotate(PI, axis=UP)
        moving_triangle.target.rotate(ang-PI/2)
        moving_triangle.target.shift(VGroup(grav_line, vert_comp, horz_comp).get_center()-moving_triangle.target.get_center())
        moving_triangle.target.set_height(vert_comp.get_length())
        moving_triangle.target.set_width(horz_comp.get_length())

        self.remove_foreground_mobject(grav_line)
        self.play(MoveToTarget(moving_triangle), run_time=2.5)
        self.wait(0.5)
        self.play(
            ang_line_label.scale, 0.5,
            ang_line_label.move_to, grav_line.get_start(),
            ang_line_label.shift, 0.45 * np.array([math.cos(-PI/2+ang/2), math.sin(-PI/2+ang/2), 0])
        )
        self.play(
            *[FadeOut(thing) for thing in 
                (*[mob for mob in moving_triangle.submobjects if mob != ang_line_arc_copy], ang_line, ang_line_arc, right_angle_plane, similar_angle_plane, similar_angle_grav)],
            self.camera_frame.move_to, VGroup(grav_line, vert_comp, horz_comp).get_center(),
            self.camera_frame.set_width, self.camera_frame.get_width()*0.6,
            FadeIn(block.FORCES["gravity"])
        )
        self.remove(grav_line)

        self.add_sound("skabretti11.m4a")
        self.wait(0.5)

        grav_eq = TexMobject("m", "g").scale(0.4).move_to(grav_line.get_center()).shift(0.3*RIGHT)
        vert_eq = TexMobject("m", "g", "\\cos", "{\\theta}").scale(0.4).rotate(PI/2).next_to(vert_comp, LEFT, buff=0.08)
        horz_eq = TexMobject("m", "g", "\\sin", "{\\theta}").scale(0.4).next_to(horz_comp, DOWN, buff=0.08)
        self.play(Write(grav_eq))
        self.wait(4)
        self.play(AnimationGroup(
            ReplacementTransform(ang_line_label.copy(), horz_eq[-1]),
            FadeIn(horz_eq[-2]), 
            lag_ratio=0.8
        ))
        self.wait(2.5)
        self.play(
            ReplacementTransform(grav_eq[0].copy(), horz_eq[0]), 
            ReplacementTransform(grav_eq[1].copy(), horz_eq[1])
        )
        self.wait(2)
        self.play(
            ReplacementTransform(ang_line_label.copy(), vert_eq[-1]),
            FadeIn(vert_eq[-2]), 
            ReplacementTransform(grav_eq[0].copy(), vert_eq[0]), 
            ReplacementTransform(grav_eq[1].copy(), vert_eq[1])
        )
        # self.play(
        #     ReplacementTransform(grav_eq[0].copy(), vert_eq[0]),
        #     ReplacementTransform(grav_eq[0].copy(), horz_eq[0]),
        #     ReplacementTransform(grav_eq[1].copy(), vert_eq[1]),
        #     ReplacementTransform(grav_eq[1].copy(), horz_eq[1]),
        #     ReplacementTransform(ang_line_label.copy(), vert_eq[-1]),
        #     ReplacementTransform(ang_line_label.copy(), horz_eq[-1]), 
        #     FadeInFrom(vert_eq[-2], 0.2*DOWN),
        #     FadeInFrom(horz_eq[-2], 0.2*LEFT),
        #     run_time=2      
        # )
        self.wait(2)

        self.play(
            ApplyMethod(self.camera_frame.restore),
            *[FadeOut(thing) for thing in (grav_eq, right_angle_grav, ang_line_arc_copy, ang_line_label)],
            *[FadeIn(thing) for thing in removed_in_grav_section],
            run_time=1.8
        )

        self.add_sound("skabretti12.m4a")
        self.wait(6.2)

        self.play(
            new_force_eq[0].shift, 1.5*DOWN,
            vec_arrows["gravity"].shift, 1.5*DOWN,
            block.FORCES["normal"].set_opacity, 0.25,
            block.FORCES["friction"].set_opacity, 0.25,
        )
        equals = TexMobject("=").next_to(new_force_eq[0], RIGHT)

        grav_vector = Matrix(["mg\\sin{\\theta}", "-mg\\cos{\\theta}"], brackets=("(", ")"), bracket_v_buff=SMALL_BUFF, v_buff=0.9)
        grav_vector.elements[0].shift((grav_vector.elements[0].get_x()-grav_vector.get_x())*LEFT)
        grav_vector.scale(0.8).set_color(self.forces["gravity"]["color"])

        norm_vector = Matrix(["0", "\\text{\\textit{Þ}}"], brackets=("(", ")"), bracket_v_buff=SMALL_BUFF, v_buff=0.9)
        norm_vector.scale(0.8).set_color(self.forces["normal"]["color"])

        fric_vector = Matrix(["-\\mu\\text{\\textit{Þ}}", "0"], brackets=("(", ")"), bracket_v_buff=SMALL_BUFF, v_buff=0.9)
        fric_vector.scale(0.8).set_color(self.forces["friction"]["color"])
        fric_vector.elements[1].shift((fric_vector.elements[1].get_x()-fric_vector.get_x())*LEFT)

        accel_vector = Matrix(["a", "0"], brackets=("(", ")"), bracket_v_buff=SMALL_BUFF, v_buff=0.9)
        accel_vector.scale(0.8)

        y_pos = (grav_vector.elements[0].get_y(), grav_vector.elements[1].get_y())
        for v in (grav_vector, norm_vector, fric_vector, accel_vector):
            v.elements[0].set_y(y_pos[0])
            v.elements[1].set_y(y_pos[1])

        grav_vector.next_to(equals, RIGHT)
        grav_vector.elements.set_opacity(0)

        horz_eq.generate_target()
        horz_eq.target.scale(1/0.4 * 0.8).move_to(grav_vector.elements[0]).set_color(self.forces["gravity"]["color"])

        vert_eq.generate_target()
        vert_eq.target.rotate(-PI/2).scale(1/0.4 * 0.8).move_to(grav_vector.elements[1], aligned_edge=RIGHT).set_color(self.forces["gravity"]["color"])

        self.play(Write(equals), Write(grav_vector))
        self.wait(1)
        self.play(
            MoveToTarget(horz_eq)
        )
        self.wait(5.5)
        self.play(
            MoveToTarget(vert_eq)
        )
        self.play(FadeOut(vert_eq), ApplyMethod(grav_vector.elements[1].set_opacity, 1))
        self.remove(horz_eq)
        grav_vector.elements[0].set_opacity(1)
        self.wait(5)

        self.play(
            ApplyMethod(grav_vector.next_to, new_force_eq[1], LEFT, buff=SMALL_BUFF),
            *[FadeOut(thing) for thing in (new_force_eq[0], vec_arrows["gravity"], equals, comps)],
            ApplyMethod(block.FORCES["normal"].set_opacity, 1),
            ApplyMethod(block.FORCES["friction"].set_opacity, 1),
        )
        self.wait(1)

        self.add_sound("skabretti13.m4a")
        self.wait(1)

        self.play(
            new_force_eq[2].shift, 1.5*DOWN,
            vec_arrows["normal"].shift, 1.5*DOWN,
            block.FORCES["gravity"].set_opacity, 0.25,
            block.FORCES["friction"].set_opacity, 0.25,
        )
        equals.next_to(new_force_eq[2])
        norm_vector.next_to(equals, RIGHT)
        self.play(Write(equals), Write(norm_vector))
        self.wait(6)
        norm_width_change = norm_vector.get_width() - new_force_eq[2].get_width() + SMALL_BUFF
        self.play(
            ApplyMethod(norm_vector.next_to, new_force_eq[3], LEFT, buff=SMALL_BUFF),
            ApplyMethod(new_force_eq[1].shift, norm_width_change*LEFT),
            ApplyMethod(grav_vector.shift, norm_width_change*LEFT),
            *[FadeOut(thing) for thing in (new_force_eq[2], vec_arrows["normal"], equals)],
            ApplyMethod(block.FORCES["gravity"].set_opacity, 1),
            ApplyMethod(block.FORCES["friction"].set_opacity, 1),
        )
        self.wait(1)

        self.add_sound("skabretti14.m4a")
        self.wait(1)

        self.play(
            new_force_eq[4].shift, 1.5*DOWN,
            vec_arrows["friction"].shift, 1.5*DOWN,
            block.FORCES["gravity"].set_opacity, 0.25,
            block.FORCES["normal"].set_opacity, 0.25,
        )
        equals.next_to(new_force_eq[4])
        fric_vector.next_to(equals, RIGHT)
        self.play(Write(equals), Write(fric_vector))
        self.wait(14)
        fric_width_change = fric_vector.get_width() - new_force_eq[4].get_width() + SMALL_BUFF
        self.play(
            ApplyMethod(fric_vector.next_to, newt_eq[1], LEFT, buff=SMALL_BUFF),
            ApplyMethod(new_force_eq[1].shift, fric_width_change*LEFT),
            ApplyMethod(new_force_eq[3].shift, fric_width_change*LEFT),
            ApplyMethod(grav_vector.shift, fric_width_change*LEFT),
            ApplyMethod(norm_vector.shift, fric_width_change*LEFT),
            *[FadeOut(thing) for thing in (new_force_eq[4], vec_arrows["friction"], equals)],
            ApplyMethod(block.FORCES["gravity"].set_opacity, 1),
            ApplyMethod(block.FORCES["normal"].set_opacity, 1),
        )
        self.wait(1)

        self.add_sound("skabretti15.m4a")
        self.wait(1)

        accel_vector_on_block = block.FORCES["total"].copy().scale(self.m, about_point=block.FORCES["total"].get_start()).set_color(WHITE)
        self.play(
            ShowCreation(accel_vector_on_block),
            newt_eq[3].shift, 1.5*DOWN,
            vec_arrows["a"].shift, 1.5*DOWN,
            block.FORCES["gravity"].set_opacity, 0.25,
            block.FORCES["normal"].set_opacity, 0.25,
            block.FORCES["friction"].set_opacity, 0.25,
        )
        self.wait(4)
        accel_vector.move_to(newt_eq[3], aligned_edge=RIGHT)
        equals.next_to(accel_vector, LEFT)
        self.play(
            ApplyMethod(VGroup(newt_eq[3], vec_arrows["a"]).next_to, equals, LEFT),
        )
        self.play(
            Write(equals), 
            Write(accel_vector)
        )
        self.wait(5)
        accel_width_change = accel_vector.get_width() - newt_eq[3].get_width() + SMALL_BUFF
        self.play(
            ApplyMethod(accel_vector.shift, 1.5*UP),
            *[ApplyMethod(thing.shift, accel_width_change*LEFT) for thing in (new_force_eq[1:4:2], grav_vector, norm_vector, fric_vector, newt_eq[1:-1])],
            *[FadeOut(thing) for thing in (newt_eq[3], vec_arrows["a"], equals)],
            ApplyMethod(block.FORCES["gravity"].set_opacity, 1),
            ApplyMethod(block.FORCES["normal"].set_opacity, 1),
            ApplyMethod(block.FORCES["friction"].set_opacity, 1),
            Uncreate(accel_vector_on_block)
        )
        self.wait(4)

        self.camera_frame.save_state()
        force_eq = VGroup(new_force_eq[1::2], grav_vector, norm_vector, fric_vector, newt_eq[1:-1], accel_vector)
        removed_before_find_acceleration = [plane, block, block.label, *[block.FORCES[k] for k in ("gravity", "normal", "friction")]]

        self.add_sound("skabretti16.m4a")
        self.wait(1)

        self.play(
            ApplyMethod(self.camera_frame.move_to, force_eq.get_center()),
            *[FadeOut(thing) for thing in removed_before_find_acceleration],
            run_time=2
        )
        self.add(new_force_eq[1::2], newt_eq[1:-1])
        self.wait(4.5)

        self.camera_frame.center()
        force_eq.center()
        
        left_plus = new_force_eq[1]
        right_plus = new_force_eq[3]
        equals = newt_eq[1]
        m_sign = newt_eq[2]

        brace = Brace(force_eq, LEFT)

        for sign in (left_plus, right_plus, equals, m_sign):
            sign.generate_target()
            sign.target = VGroup(
                sign.copy().scale(0.8).set_y(y_pos[0]),
                sign.copy().scale(0.8).set_y(y_pos[1])
            ).set_x(sign.get_x())

        self.add_sound("skabretti17.m4a")
        self.wait(1)

        self.play(
            *[FadeOut(brackets) for brackets in [vec.brackets for vec in (grav_vector, norm_vector, fric_vector, accel_vector)]],
            *[MoveToTarget(sign) for sign in (left_plus, right_plus, equals, m_sign)],
            GrowFromCenter(brace),
            run_time=2.5
        )
        self.wait(0.6)

        fric_tex = TexMobject("-", "\\mu\\text{\\textit{Þ}}").scale(0.8).set_color(self.forces["friction"]["color"]) \
            .move_to(fric_vector.elements[0].get_center())
        self.remove(fric_vector.elements[0])
        self.add(fric_tex)

        x_eq = TexMobject("mg", "\\sin{\\theta}", "-", "\\mu", "\\text{\\textit{Þ}}", "=", "m", "a").scale(0.8)
        y_eq = TexMobject("-", "mg", "\\cos{\\theta}", "+", "\\text{\\textit{Þ}}", "=", "0").scale(0.8)
        x_eq.set_y(y_pos[0])
        y_eq.set_y(y_pos[1])

        for g in (x_eq[0], x_eq[1], y_eq[0], y_eq[1], y_eq[2]):
            g.set_color(self.forces["gravity"]["color"])
        y_eq[4].set_color(self.forces["normal"]["color"])
        x_eq[3].set_color(self.forces["friction"]["color"])
        x_eq[4].set_color(self.forces["friction"]["color"])

        self.play(AnimationGroup(
            AnimationGroup(
                *[FadeOut(thing) for thing in (left_plus[0], right_plus, m_sign[1], norm_vector.elements[0], fric_vector.elements[1])],
                lag_ratio=0.1
            ),
            AnimationGroup(
                *[ReplacementTransform(old, new) for old, new in zip(
                    [fric_tex[0], left_plus[1], fric_tex[1], norm_vector.elements[1], 
                        equals[0], equals[1], m_sign[0], accel_vector.elements[1], accel_vector.elements[0]], 
                    [x_eq[2], y_eq[3], VGroup(x_eq[3], x_eq[4]), y_eq[4], 
                        x_eq[5], y_eq[5], x_eq[6], y_eq[6], x_eq[7]]
                )],
                *[ApplyMethod(old.move_to, new.get_center()) for old, new in zip(grav_vector.elements, [VGroup(x_eq[0], x_eq[1]), VGroup(y_eq[0], y_eq[1], y_eq[2])])],
                ApplyMethod(brace.next_to, VGroup(x_eq, y_eq), LEFT),
                run_time=1.5
            ),
            lag_ratio=0.7,
        ))
        self.remove(*grav_vector.elements)
        self.add(x_eq, y_eq)

        self.play(AnimationGroup(
            AnimationGroup(WiggleOutThenIn(x_eq[-4]), WiggleOutThenIn(y_eq[-3])),
            WiggleOutThenIn(x_eq[-1]),
            lag_ratio=0.75
        ))

        self.add_sound("skabretti18.m4a")
        self.wait(7)
        
        self.play(
            self.camera_frame.set_width, y_eq.get_width()*3,
            self.camera_frame.move_to, y_eq.get_center(),
            *[ApplyMethod(thing.set_opacity, 0.25) for thing in (x_eq, brace)],
        )
        self.wait(7)

        self.add_sound("skabretti19.m4a")
        self.wait(1)

        self.play(AnimationGroup(
            *[FadeOut(y_eq[i]) for i in (0, 3, 6)],
            MoveAlongPath(VGroup(y_eq[1], y_eq[2]), ArcBetweenPoints(VGroup(y_eq[1], y_eq[2]).get_center(), y_eq[6].get_left()+VGroup(y_eq[1], y_eq[2]).get_width()/2*RIGHT)),
            lag_ratio=0.2,
            run_time=1.8
        ))
        y_eq = VGroup(y_eq[4], y_eq[5], y_eq[1], y_eq[2])
        y_eq.generate_target()
        y_eq.target.align_to(x_eq, LEFT)
        y_eq.target[-2:].set_color(self.forces["normal"]["color"])

        self.play(
            MoveToTarget(y_eq)
        )
        self.wait(4)
        self.play(
            self.camera_frame.set_y, VGroup(x_eq, y_eq, brace).get_y(),
            *[ApplyMethod(thing.set_opacity, 1) for thing in (x_eq, brace)]
        )
        self.wait(0.5)

        mu_buff = 0.05
        x_eq_width_change = y_eq[-2:].get_width() - x_eq.get_parts_by_tex("Þ").get_width() + mu_buff
        self.play(
            ApplyMethod(x_eq[-3:].shift, x_eq_width_change*RIGHT),
            ApplyMethod(y_eq[-2:].move_to, x_eq[-4].copy().shift(mu_buff*RIGHT), {"aligned_edge": UL}),
            *[FadeOut(thing) for thing in (y_eq[:-2], brace, x_eq.get_parts_by_tex("Þ"))],
            self.camera_frame.move_to, x_eq.get_center()+x_eq_width_change/2*RIGHT,
            run_time=1.5
        )
        self.wait(1)

        mg_out_width_change = y_eq[-2].get_right()[0]-x_eq[3].get_right()[0]

        paren_buff = 0.075
        left_paren = TexMobject("(").scale(0.8).next_to(x_eq[1], LEFT, buff=paren_buff).shift(mg_out_width_change*RIGHT)
        right_paren = TexMobject(")").scale(0.8).next_to(y_eq[-1], RIGHT, buff=paren_buff)

        self.add_foreground_mobject(x_eq[0])
        x_eq[0].generate_target()
        x_eq[0].target.shift((left_paren.get_width()+paren_buff-mg_out_width_change)*LEFT)
        y_eq[-2].generate_target()
        y_eq[-2].target.move_to(x_eq[0].target).set_opacity(0.25).set_color(self.forces["gravity"]["color"])

        self.play(AnimationGroup(
            AnimationGroup(
                ApplyMethod(x_eq[:4].shift, (mg_out_width_change)*RIGHT),
                ApplyMethod(x_eq[-3:].shift, (right_paren.get_width()+paren_buff)*RIGHT),
                MoveToTarget(x_eq[0]),
                ClockwiseTransform(y_eq[-2], y_eq[-2].target),
                ApplyMethod(
                    self.camera_frame.shift, 
                    ((x_eq[0].target.get_center()[0]-x_eq[0].get_center()[0]) - (right_paren.get_width()+paren_buff))*RIGHT
                ),
                run_time=2
            ),
            AnimationGroup(
                FadeIn(left_paren),
                FadeIn(right_paren)
            ),
            lag_ratio=0.5,
        ))
        self.play(FadeOut(y_eq[-2]), run_time=0.5)
        self.wait(1)

        right_mg = TexMobject("m", "g").scale(0.8).match_style(x_eq[0]).move_to(x_eq[0])
        self.add(right_mg)
        self.remove(x_eq[0])
        self.play(
            FadeOutAndShift(right_mg[0], 0.3*DOWN),
            FadeOutAndShift(x_eq[-2], 0.3*DOWN),
            ApplyMethod(x_eq[-1].shift, (x_eq[-2].get_width()+0.01)*LEFT)
        )
        self.wait(0.75)

        for s in x_eq[-3], x_eq[-1]:
            s.generate_target()
            s.target.set_x(VGroup(right_mg[1], x_eq[1:]).get_x()-s.get_x())

        self.play(
            CounterclockwiseTransform(x_eq[-3], x_eq[-3].target),
            CounterclockwiseTransform(x_eq[-1], x_eq[-1].target),
            ApplyMethod(self.camera_frame.set_x, VGroup(x_eq[-1].target, right_paren).get_x()),
            run_time=1.25
        )

        self.wait(10)

        plane = InclinedPlane(
            angle=self.theta, 
            show_angle_arc=False, 
            show_angle_label=False,
            arc_radius=1.5
        ).to_corner(DL)

        block = MovingBlock(
            mass=self.m, 
            start_point=plane.get_vertices()[1], 
            end_point=plane.get_vertices()[2], 
            displacement=1,
            label_scale=1.25,
            unit_height=0.9
        )

        for f in self.forces.keys():
            self.forces[f]["angle"] -= ang

        self.make_forces(block)
        accel_vector_on_block = block.FORCES["total"].copy().scale(self.m, about_point=block.FORCES["total"].get_start()).set_color(WHITE)
        accel_vector_on_block.add_updater(lambda a: a.move_to(block.get_critical_point_on_block(RIGHT), aligned_edge=UL))

        acc_eq = TexMobject("a", "=", "g", "(", "\\sin{\\theta}", "-", "\\mu", "\\cos{\\theta}", ")").scale(0.8)
        for i, j in zip(range(len(acc_eq)), [x_eq[-1], x_eq[-3], right_mg[-1], left_paren, x_eq[1], x_eq[2], x_eq[3], y_eq[3], right_paren]):
            acc_eq[i].move_to(j).match_style(j)
        self.add(acc_eq)
        self.remove(*x_eq.submobjects, *y_eq.submobjects, left_paren, right_paren, *right_mg.submobjects)

        acc_eq_cam_offset = acc_eq.get_center() - self.camera_frame.get_center()
        acc_eq.next_to(accel_vector_on_block)
        acc_eq_block_offset = acc_eq.get_center() - block.get_center()
        acc_eq.add_updater(lambda eq: eq.move_to(block.get_center() + acc_eq_block_offset))
        
        self.camera_frame.move_to(acc_eq.get_center() - acc_eq_cam_offset)
        self.camera_frame.generate_target()
        self.camera_frame.target.move_to(ORIGIN).set_width(FRAME_WIDTH)

        self.play(MoveToTarget(self.camera_frame), FadeIn(plane), FadeIn(block), FadeIn(accel_vector_on_block))
        self.wait(1)

        acc_eq.add_updater(lambda eq: eq.next_to(accel_vector_on_block))

        block.start_move(a)
        self.wait(slide_time/2)
        vel = block.velocity
        acc_eq.clear_updaters()
        self.play(FadeOutAndShift(acc_eq, vel*np.array([math.cos(-ang), math.sin(-ang), 0])), rate_func=linear)
        self.wait(slide_time/2)



    def make_forces(self, block):
        for f in self.forces.keys():
            block.make_force_vector(
                f, 
                self.forces[f]["size"] * np.array([math.cos(self.forces[f]["angle"]), math.sin(self.forces[f]["angle"]), 0]), 
                self.forces[f]["origin"], 
                color=self.forces[f]["color"],
            )
            block.make_force_label(
                f,
                self.forces[f]["label_text"], 
                self.forces[f]["label_direction"],
                scale=min(0.8, self.forces[f]["size"]),
                color=self.forces[f]["color"]
            )

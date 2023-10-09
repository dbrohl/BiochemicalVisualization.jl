struct CubicB

end

function (spline::CubicB)()

end

function from_protein_chain(chain)
    n_ribbons = 9
    ribbon_width = Dict(SecondaryStructure.HELIX => 3, SecondaryStructure.SHEET => 3, SecondaryStructure.NONE => 1)
    guide_points_per_segment = 10

    # template for debug output
    colors = [(255, 0, 0), (255, 255, 0), (0, 255, 0), (0, 255, 255), (0, 0, 255), (255, 0, 255)]
    sphere_radius = 0.2
    sphere_mesh = simplexify(Sphere{3, Float64}((0,0,0), sphere_radius))
    debug_mesh = ColoredMesh(sphere_mesh, (0, 0, 0))

    c_alphas = filter(x -> x.element==Elements.C && x.name=="CA", atoms(chain))
    oxygens = filter(x -> x.element==Elements.O && x.name=="O", atoms(chain))
    @assert length(c_alphas)==length(oxygens)
    structures = [SecondaryStructure.NONE for a in c_alphas] # TODO 

    for a in c_alphas
        debug_mesh = BiochemicalVisualization.merge(debug_mesh, ColoredMesh(Translate(Float64.(a.r)...)((sphere_mesh)), (100, 100, 100)))
    end



    guide_points = Array{Float64, 3}(undef, guide_points_per_segment, 3, length(c_alphas)-1)
    for i=1:length(c_alphas)-1
        A = c_alphas[i+1].r - c_alphas[i].r
        B = oxygens[i].r - c_alphas[i].r
        C = cross(A, B)
        D = cross(C, A)

        C = C/norm(C)
        D = D/norm(D)

        P = c_alphas[i].r + 0.5 * A
        if(structures[i]==SecondaryStructure.HELIX || structures[i+1]==SecondaryStructure.HELIX) # TODO correct condition?
            # TODO translate
        end
        D *= 0.5 * ribbon_width[structures[i]]

        for (j, t) in enumerate(range(-1, 1, guide_points_per_segment))
            guide_points[j, :, i] = P + t*D
            
            m = Translate(guide_points[j, :, i]...)((sphere_mesh))
            m = ColoredMesh(m, colors[mod1(i, length(colors))])
            debug_mesh = BiochemicalVisualization.merge(debug_mesh, m)
        end


    end

    # a = 1
    # while a+3<=size(guide_points, 3)
    #     debug_mesh = BiochemicalVisualization.merge(debug_mesh, b_spline_quadrupel(guide_points[1, :, a], guide_points[1, :, a+1], guide_points[1, :, a+2], guide_points[1, :, a+3]))
    #     a+=1
    # end
    debug_mesh = BiochemicalVisualization.merge(debug_mesh, b_spline_quadrupel(guide_points[1, :, 1], guide_points[1, :, 2], guide_points[1, :, 3], guide_points[1, :, 4]))

    export_mesh_to_ply("ribbon_debug.ply", debug_mesh)
end

# function b_spline_quadrupel(P1, P2, P3, P4)
#     debug_mesh = ColoredMesh(SimpleMesh([(0, 0, 0), (1, 1, 1)], [connect((1, 2))]), (255, 0, 0)) 

#     N = 10
#     S = [6/N^3 0     0   0
#          6/N^3 2/N^2 0   0
#          1/N^3 1/N^2 1/N 0
#          0     0     0   1]
#     B = 1/6 * [-1 3 -3 1
#                 3 -6 3 0
#                 -3 0 3 0
#                 1 4 1 0]
#     G = [P1 P2 P3 P4]'
#     G = [G [1
#     1
#     1
#     1]]
#     log_info(misc, "P1", P1)
#     log_info(misc, "P2", P2)
#     log_info(misc, "P3", P3)
#     log_info(misc, "P4", P4)
#     log_info(misc, "G", G)
#     M = S * B * G

#     log_info(misc, "M", M)

#     current_coords = M[4, 1:3]/M[4,4]
#     for k=1:N
#         for i=4:2
#             for j =1:4
#                 M[i,j] += M[i-1, j]
#             end
#         end
#         next_coords = M[4, 1:3]/M[4,4]
#         log_info(misc, k, next_coords)
#         points = [Tuple(current_coords), Tuple(next_coords)]
#         debug_mesh = BiochemicalVisualization.merge(debug_mesh, ColoredMesh(SimpleMesh(points, [connect((1, 2))]), (255, 0, 0)))

#         current_coords = next_coords
#     end
#     return debug_mesh
# end





# function Knot(i)
#     if(i < knotK)
#         Knot = 0
#     elseif(i > knotN)
#         Knot = knotN — knotK + 2
#     else
#         Knot = i — knotK + 1
#     end
#     return Knot
# end

# function NBlend(i, k, u)
#     if(k==1)
#         v=0
#         if (Knot(i)<=u && u<Knot(i + 1))
#             v = 1
#         end
#     else 
#         v=0
#         t = Knot(i + k — 1) — Knot(i)
#         if(t!=0)
#             v = (u - Knot(i)) * NBlend(i, k - 1, u) / t
#         end
#         t = Knot(i + k) — Knot(i + 1)
#         if(t!=0)
#             v += (Knot(i + k) — u)* NBlend(i + 1, k — 1, u) / t
#         end
#     end
#     return v

# end;

# function BSpline(x, y, z, u, n, k, p )
#     knotK = k
#     knotN = n
    
#     x=0
#     y=0
#     z=0
#     for i=1:n
#         b = NBlend(i, k, u)
#         x += p[i, 1]*b
#         y += p[i, 2]*b
#         z += p[i, 3]*b
#     end
# end
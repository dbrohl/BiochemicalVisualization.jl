# create reduced surface
    # create list of spheres(atom positions and radii+probe_radius (+ ggf. offset))
    # preprocessing
        # determine bounding box of atom centers (without vdw/probe radius)
        # determine max radius
        # create hash grid (Sanner et al use a tree instead of a grid. )

#TODO separate methode für chains? sollten chains als verschiedene Oberflächen behandelt werden?


function create_reduced_surface(
    ac::System{T}, probe_radius::T) where {T<:Real}
    if natoms(ac)==0
        return # TODO better error
    end

    epsilon = 1e-4;

    first_atom = first(eachatom(ac))
    max_radius = 0
    bounding_box = [first_atom.r;; first_atom.r] # [x/y/z, min/max]
    println(bounding_box)

    for atom in eachatom(ac)
        max_radius = max(max_radius, atom.radius==0 ? 1 : atom.radius)

        bounding_box[:, 1] = min.(bounding_box[:, 1], atom.r)
        bounding_box[:, 2] = max.(bounding_box[:, 2], atom.r)
    end

    # dist = 2*(max_radius + probe_radius)
    # n_boxes = (bounding_box[:, 2] - bounding_box[:, 1])/dist + 5 #TODO warum +5?

    # hashgrid = CellList(map(atom -> atom.r, eachatom(ac)), dist)

    # HashGridBox3<Position>* box;
    # HashGridBox3<Position>::ConstBoxIterator b;
    # HashGridBox3<Position>::ConstDataIterator d;

    # std::list<Position> to_delete;
    # Size num_deleted = 0;
    # Vector3 pos;
    # for (Position i = 0; i < rs_->number_of_atoms_; i++)
    # {
    #     pos.set(rs_->atom_[i].p.x, rs_->atom_[i].p.y, rs_->atom_[i].p.z);

    #     // remove atoms that are fully contained in another atom
    #     double radius_i = rs_->atom_[i].radius;
    #     bool too_close = false;
    #     box = grid.getBox(pos);
    #     for (b = box->beginBox(); b != box->endBox() && !too_close; b++)
    #     {
    #         for (d = b->beginData(); d != b->endData() && !too_close; d++)
    #         {
    #             double radius_d = rs_->atom_[*d].radius;

    #             // our algorithm becomes instable somewhere if two atoms are too close...
    #             // TODO: fix it so this safe guard is no longer necessary
    #             if (Maths::isLessOrEqual(rs_->atom_[i].p.getDistance(rs_->atom_[*d].p), 0.05*std::max(radius_d, radius_i)))
    #             {
    #                 too_close = true;
    #                 to_delete.push_back(i);
    #                 num_deleted++;
    #                 if (radius_i > radius_d)
    #                 {
    #                     rs_->atom_[*d].p = rs_->atom_[i].p;
    #                     rs_->atom_[*d].radius = rs_->atom_[i].radius;
    #                 }
    #             }
    #         }
    #     }
        
    #     if (!too_close)
    #         grid.insert(pos, i-num_deleted);
    # }

    # for atom in eachatom(ac)
    #     for neighbor in hashgrid.get(atom)

    #         if (distance(atom.r, neighbor.r) - 0.05*max(neighbor.radius, atom.radius)) < epsilon

    #             too_close = true
                
    #         if too_close
    #             break
    #         end

    #     end
    # end

end
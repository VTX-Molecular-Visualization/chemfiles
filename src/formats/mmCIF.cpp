// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/parse.hpp"
#include "chemfiles/string_view.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/utils.hpp"

#include "chemfiles/Atom.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/FormatMetadata.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"

// WARNING UGLY HACK!
#include "chemfiles/formats/PDB.hpp"

#include "chemfiles/formats/mmCIF.hpp"

using namespace chemfiles;

template<> const FormatMetadata& chemfiles::format_metadata<mmCIFFormat>() {
    static FormatMetadata metadata;
    metadata.name = "mmCIF";
    metadata.extension = ".mmcif";
    metadata.description = "Crystallographic Information Framework files for MacroMolecules";
    metadata.reference = "http://mmcif.wwpdb.org/";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = true;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = true;
    metadata.residues = true;
    return metadata;
}

/// CIF files store which digits are insignificant, we need to remove this
static double cif_to_double(std::string line);

void mmCIFFormat::init_() {
    if (file_.mode() == File::WRITE) {
        return;
    }

    if (file_.mode() == File::APPEND) {
        throw file_error("cannot open mmCIF files in append ('a') mode");
    }

    Vector3D lengths;
    Vector3D angles = {90, 90, 90};

    reader_meta_data_ = ReaderMetaData();

    bool in_loop = false;
    
    while (!file_.eof())
    {
        auto line = file_.readline();
        reader_meta_data_.line = line;

        if (line.find("loop_") != std::string::npos) {
            in_loop = true;
            continue;
        }

        if (line.empty() || line[0] == '#') {
            in_loop = false;
            continue;
        }

        auto line_split = split(line, ' ');
        if (line_split.size() > 1 && line[0] == '_') 
        {
            in_loop = false;
        }

        if (!in_loop )
        { 
            if (line_split[0] == "_cell_length_a" || line_split[0] == "_cell.length_a") {
                lengths[0] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_cell_length_b" || line_split[0] == "_cell.length_b") {
                lengths[1] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_cell_length_c" || line_split[0] == "_cell.length_c") {
                lengths[2] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_cell_angle_alpha" || line_split[0] == "_cell.angle_alpha") {
                angles[0] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_cell_angle_beta" || line_split[0] == "_cell.angle_beta") {
                angles[1] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_cell_angle_gamma" || line_split[0] == "_cell.angle_gamma") {
                angles[2] = cif_to_double(line_split[1].to_string());
            }
            if (line_split[0] == "_entry.id") {
                pdb_idcode_ = line_split[1].to_string();
            }

            if (line_split[0] == "_struct.title") 
            {
                if ( line_split.size() > 1 )
                {
                    auto struct_name = line.substr(13);
                    struct_name = trim(struct_name);

                    if ( struct_name[0] == '\'' )
                        struct_name = struct_name.substr(1, struct_name.size() - 2);
                   
                    name_ = struct_name.to_string();
                }
                else
                {
                    line = file_.readline();
                    trim(line);

                    if ( line[0] == ';' )  //Multi-line title
                    {
                        name_ = read_multi_line();
                    }
                    else if ( line[0] == '\'' )    //Single line in quote
                    {
                        name_ = line.substr(1, line.size() - 2).to_string();
                    }
                    else //Single line
                    {
                        name_ = line.to_string();
                    }
                }
            }
        }
        else // if (in_loop)
        { 
            if ( line_split[0].find("_atom_site.") != std::string::npos ) {
                init_atom_site();
            }
        }

        // Can be in loop or not
        if ( line_split[0].find("_pdbx_struct_assembly_gen.") != std::string::npos ) {
            addCategory("_pdbx_struct_assembly_gen", in_loop);
        }
        
        if ( line_split[0].find("_pdbx_struct_oper_list.") != std::string::npos ) {
            addCategory("_pdbx_struct_oper_list", in_loop);
        }
    }

    if ( category_map_.find("_atom_site") == category_map_.end() ) {
        throw format_error("could not find _atom_site loop in '{}'", file_.path());
    }

    cell_ = UnitCell(lengths, angles);
}

void mmCIFFormat::init_atom_site()
{
    auto & category = category_map_["_atom_site"];
    category.is_loop = true;

    fill_loop_properties("_atom_site", category.property_map);

    auto atom_site_map = category.property_map;
    category.position_data_start = reader_meta_data_.position;

    // After this block ends, we have the start of coordinates
    steps_positions_.push_back(reader_meta_data_.position);

    if ( atom_site_map.find("type_symbol") == atom_site_map.end() ) {
        throw format_error("could not find _atom_site.type_symbol in '{}'", file_.path());
    }
    if ( atom_site_map.find("Cartn_x") == atom_site_map.end() ) {
        throw format_error("could not find _atom_site.Cartn_x in '{}'", file_.path());
    }

    // Do we have a special extension for multiple modes?
    auto model_position = atom_site_map.find("pdbx_PDB_model_num");
    if ( model_position == atom_site_map.end() ) {
        // If not, we are done
        file_.seekpos(steps_positions_[0]);
        return;
    }

    auto line = file_.readline();
    reader_meta_data_.line = line;

    // Ok, let's look at the sites now to note where models start
    auto last_position = parse<size_t>(split(line, ' ')[model_position->second]);

    do
    {
        reader_meta_data_.position = file_.tellpos();
        line = file_.readline();
        //reader_meta_data_.line = line;

        // a break in the text ends the models
        if ( line.empty() || line == "line_" || line[0] == '#' ) {
            break;
        }

        auto line_split = split(line, ' ');
        size_t current_position = parse<size_t>(line_split[model_position->second]);

        if ( current_position != last_position ) {
            steps_positions_.push_back(reader_meta_data_.position);
            last_position = current_position;
        }
    } while ( !file_.eof() );

    category.position_data_end = reader_meta_data_.position;
}
void mmCIFFormat::addCategory(const std::string& category_name, bool is_loop)
{
    auto & category = category_map_[category_name];
    category.is_loop = is_loop;

    if ( is_loop )
    {
        fill_loop_properties(category_name, category.property_map);
       
        category.position_data_start = reader_meta_data_.position;

        while ( !file_.eof() )
        {
            reader_meta_data_.position = file_.tellpos();
            auto line = file_.readline();
            if ( line.empty() || line[0] == '#' )
            {
                break;
            }
        }

        category.position_data_end = reader_meta_data_.position;
    }
    else
    {
        category.position_data_start = file_.tellpos() - (reader_meta_data_.line.size() + 1);

        const size_t category_name_length = category_name.size() + 1;    //Size + separator
        size_t property_index = 0;

        // Fill attributes
        auto line = reader_meta_data_.line;

        do
        {
            auto line_split = split(line, ' ');

            if ( line.find(category_name) != std::string::npos )
            {
                auto label = trim(line_split[0]).substr(category_name_length).to_string();
                category.property_map[label] = property_index++;

                reader_meta_data_.position = file_.tellpos();
                line = file_.readline();
                reader_meta_data_.line = line;

                continue;
            }

            // Anything else ends the block
            break;

        } while ( !file_.eof() );

        category.position_data_end = reader_meta_data_.position + 1;
    }
}
void mmCIFFormat::fill_loop_properties(const std::string& category_name, std::map<std::string, size_t>& property_map)
{
    const size_t category_name_length = category_name.size() + 1;    //Size + separator
    size_t property_index = 0;

    // Fill attributes
    auto line = reader_meta_data_.line;

    do
    {
        auto line_split = split(line, ' ');

        if ( line.find(category_name) != std::string::npos )
        {
            auto atom_label = trim(line).substr(category_name_length).to_string();
            property_map[atom_label] = property_index++;

            reader_meta_data_.position = file_.tellpos();
            line = file_.readline();
            reader_meta_data_.line = line;

            continue;
        }

        // Anything else ends the block
        break;

    } while ( !file_.eof() );
}

void mmCIFFormat::read_item_properties(const mmCIFCategoryHeader& header, std::vector<std::string>& properties)
{
    properties.clear();

    read_property_line(properties);

    while ( properties.size() < header.property_map.size() )
    {
        if ( reader_meta_data_.line.empty() || reader_meta_data_.line == "loop_" || reader_meta_data_.line[0] == '#' ) {
            break;
        }

        read_property_line(properties);
    }
}
void mmCIFFormat::read_property_line(std::vector<std::string>& properties)
{
    reader_meta_data_.position = file_.tellpos();
    auto line = file_.readline();
    reader_meta_data_.line = line;

    size_t wordStart = 0;
    size_t current = 0;
    bool inString = false;

    while ( current < line.size() )
    {
        auto currentChar = line[current];
        if ( currentChar == '\'' )
        {
            if ( inString )
            {
                properties.emplace_back(line.substr(wordStart, current - wordStart).to_string());
                inString = false;
            }
            else
            {
                inString = true;
            }

            current++;
            wordStart = current;
        }
        else if ( currentChar == ' ' && !inString )
        {
            if ( wordStart != current )
            {
                properties.emplace_back(line.substr(wordStart, current - wordStart).to_string());
            }

            current++;
            wordStart = current;
        }
        else
        {
            current++;
        }
    }
}

size_t mmCIFFormat::nsteps() {
    return steps_positions_.size();
}
void mmCIFFormat::read_step(const size_t step, Frame& frame) {
    assert(step < steps_positions_.size());
    file_.seekpos(steps_positions_[step]);
    read(frame);
}
void mmCIFFormat::read(Frame& frame) {
    map_residues_indexes.clear();
    residues_.clear();
    frame.set_cell(cell_);

    if (!name_.empty()) {
        frame.set("name", name_);
    }

    if (!pdb_idcode_.empty()) {
        frame.set("pdb_idcode", pdb_idcode_);
    }

    read_atom_site(frame);
    for ( const auto& residue : residues_ )
    {
        frame.add_residue(residue);
    }

    auto atom_site_map = category_map_["_atom_site"].property_map;
    auto model_position = atom_site_map.find("pdbx_PDB_model_num");

    // Only link if we are reading mmCIF
    if ( model_position != atom_site_map.end() )
    {
        // Cross format talk! Forgive me!
        PDBFormat::link_standard_residue_bonds(frame);
    }

    fill_assembly();

    if (assembly_.assembly_generators.size() > 0 )
        apply_symmetry(frame, assembly_.assembly_generators.begin()->assembly_id);
}

void mmCIFFormat::read_atom_site(Frame& frame)
{   // The following map operations can be moved to the constructor
    // and stored with optional. I don't know which solution is better.

    auto atom_site = category_map_["_atom_site"];
    auto atom_site_map = atom_site.property_map;

    // These are required for atoms
    auto type_symbol = atom_site_map.at("type_symbol");

    // This has two names...
    auto label_atom_id = atom_site_map.find("label_atom_id");
    if ( label_atom_id == atom_site_map.end() )
    {
        label_atom_id = atom_site_map.find("label");
    }

    // Other atom properties
    auto group_pdb = atom_site_map.find("group_PDB");
    auto label_alt_id = atom_site_map.find("label_alt_id");
    auto formal_charge = atom_site_map.find("formal_charge");

    // Positions
    auto cartn_x = atom_site_map.at("Cartn_x");
    auto cartn_y = atom_site_map.at("Cartn_y");
    auto cartn_z = atom_site_map.at("Cartn_z");

    // Residue properties
    auto label_comp_id = atom_site_map.find("label_comp_id");
    auto label_asym_id = atom_site_map.find("label_asym_id");
    auto auth_asym_id = atom_site_map.find("auth_asym_id");
    auto label_seq_id = atom_site_map.find("label_seq_id");
    auto auth_seq_id = atom_site_map.find("auth_seq_id");
    auto label_entity_id = atom_site_map.find("label_entity_id");

    auto model_position = atom_site_map.find("pdbx_PDB_model_num");

    file_.seekpos(steps_positions_[0]);
    reader_meta_data_.position = steps_positions_[0];

    size_t last_step_position = 0;
    if ( model_position != atom_site_map.end() )
    {
        auto line = file_.readline();
        last_step_position = parse<size_t>(split(line, ' ')[model_position->second]);
        // Reset file position so that the loop below can start by reading the
        // first line
        file_.seekpos(steps_positions_[0]);
    }

    while ( !file_.eof() )
    {
        auto line = file_.readline();
        reader_meta_data_.position = file_.tellpos();

        if ( file_.tellpos() >= atom_site.position_data_end )
        {
            break;
        }

        auto line_split = split(line, ' ');
        if ( line_split.size() != atom_site_map.size() )
        {
            throw format_error("line '{}' has {} items not {}",
                               line, line_split.size(), atom_site_map.size()
            );
        }

        size_t current_step_position = 0;
        if ( model_position != atom_site_map.end() )
        {
            current_step_position = parse<size_t>(line_split[model_position->second]);
        }

        if ( current_step_position != last_step_position )
        {
            break;
        }

        Atom atom(
            line_split[label_atom_id->second].to_string(),
            line_split[type_symbol].to_string()
        );

        if ( label_alt_id != atom_site_map.end() &&
            line_split[label_alt_id->second] != "." )
        {
            atom.set("altloc", line_split[label_alt_id->second].to_string());
        }

        if ( formal_charge != atom_site_map.end() )
        {
            atom.set_charge(cif_to_double(line_split[formal_charge->second].to_string()));
        }

        auto x = cif_to_double(line_split[cartn_x].to_string());
        auto y = cif_to_double(line_split[cartn_y].to_string());
        auto z = cif_to_double(line_split[cartn_z].to_string());
        frame.add_atom(std::move(atom), Vector3D(x, y, z));


        if ( label_comp_id == atom_site_map.end() || label_asym_id == atom_site_map.end() )
        {
            continue;
        }

        auto atom_id = frame.size() - 1;
        int64_t resid = 0;

        auto residue_id_parameter = auth_seq_id != atom_site_map.end() ? auth_seq_id : label_seq_id;
        auto resid_text = line_split[residue_id_parameter->second];

        try
        {
            if ( resid_text == "." )
            { // In this case, we need to use the entity id
                resid = parse<int64_t>(line_split[label_entity_id->second]);
            }
            else
            {
                resid = parse<int64_t>(line_split[residue_id_parameter->second]);
            }
        }
        catch ( const Error& e )
        {
            throw format_error("invalid CIF residue or entity numeric: {}", e.what());
        }

        auto chainid = line_split[label_asym_id->second].to_string();

        if ( map_residues_indexes.find({ chainid, resid }) == map_residues_indexes.end() )
        {
            auto name = line_split[label_comp_id->second];
            Residue residue(name.to_string(), resid);
            residue.add_atom(atom_id);

            // This will be saved as a string on purpose to match MMTF
            if ( label_asym_id != atom_site_map.end() )
            {
                residue.set("chainid", chainid);
            }

            if ( auth_asym_id != atom_site_map.end() )
            {
                residue.set("chainname", line_split[auth_asym_id->second].to_string());
            }

            if ( group_pdb != atom_site_map.end() )
            {
                residue.set("is_standard_pdb", line_split[group_pdb->second] == "ATOM");
            }

            map_residues_indexes.emplace(std::make_pair(chainid, resid), residues_.size());
            residues_.emplace_back(std::move(residue));
        }
        else
        {
            // Just add this atom to the residue
            residues_[map_residues_indexes.at({ chainid, resid })].add_atom(atom_id);
        }
    }
}

void mmCIFFormat::apply_symmetry(Frame& frame, const std::string & assembly_id)
{
    const auto original_size = frame.size();
    const auto original_residue_size = frame.topology().residues().size();
    const auto original_bond_size = frame.topology().bonds().size();

    using bond_w_order = std::pair<Bond, Bond::BondOrder>;
    std::vector<bond_w_order> bonds_to_add;

    size_t assembly_generator_count = 0;
    for (auto assemblyGenerator : assembly_.assembly_generators )
    {
        if ( assemblyGenerator.assembly_id != assembly_id )
            continue;

        assembly_generator_count++;

        // Skip identity operation
        if ( assemblyGenerator.operations.size() == 1 )
        {
            const auto& operation = assembly_.assembly_operations[assemblyGenerator.operations[0]];
            if ( operation.rotation == Matrix3D::unit() && operation.translation == Vector3D() )
            {
                continue;
            }
        }

        std::vector<size_t> old_to_sym(original_size, 0);
        std::vector<Residue> residues_to_add;
        residues_to_add.reserve(original_residue_size);

        for ( size_t residue_index = 0 ; residue_index < original_residue_size ; ++residue_index)
        {
            const auto& residue = frame.topology().residues()[residue_index];
            auto chainID = residue.get("chainid");

            if ( !chainID || assemblyGenerator.targets.find(chainID->as_string()) == assemblyGenerator.targets.end() ) { 
                continue;
            }

            // Copy over everything except the atoms
            auto new_residue = Residue(residue.name(), *residue.id());
            for ( auto& prop : residue.properties() )
            {
                new_residue.set(prop.first, prop.second);
            }

            // Avoid using this chain in future symmetry operations
            const std::string newChainID = chainID->as_string() + "-" + std::to_string(assembly_generator_count);
            new_residue.set("chainid", newChainID);
            // Rename chain with assembly ID
            std::string newChainName = residue.get("chainname")->as_string() + "-" + std::to_string( assembly_generator_count);
            new_residue.set("chainname", newChainName);

            for ( auto atom_id : residue )
            {
                // Ensure that the current atom is not a result of a symmetry operation
                if ( atom_id >= original_size )
                {
                    continue;
                }

                auto new_atom = frame[atom_id];
                auto new_position = frame.positions()[atom_id];
               
                for ( auto operation_str : assemblyGenerator.operations )
                {
                    const auto & operation = assembly_.assembly_operations[operation_str];

                    if ( operation.rotation == Matrix3D::unit() && operation.translation == Vector3D() )
                    {
                        continue;
                    }

                    new_atom = frame[atom_id];
                    new_position = operation.rotation * new_position + operation.translation;
                }

                frame.add_atom(std::move(new_atom), std::move(new_position));
                new_residue.add_atom(frame.size() - 1);
                old_to_sym[atom_id] = frame.size() - 1;
            }

            residues_to_add.emplace_back(new_residue);
        }

        residues_to_add.shrink_to_fit();

        for ( auto&& residue : residues_to_add )
        {
            frame.add_residue(std::move(residue));
        }

        for ( size_t i = 0; i < original_bond_size; ++i )
        {
            auto& bond = frame.topology().bonds()[i];

            // bonds should be sorted so that when we hit original size, we're done
            if ( bond[0] >= original_size || bond[1] >= original_size )
            {
                break;
            }

            auto new_bond_0 = old_to_sym[bond[0]];
            auto new_bond_1 = old_to_sym[bond[1]];

            // zero simply means that the atom has no partner
            if ( new_bond_0 == 0 || new_bond_1 == 0 )
            {
                continue;
            }

            bonds_to_add.push_back({ {new_bond_0, new_bond_1}, frame.topology().bond_orders()[i] });
        }
    }

    for ( auto& bond : bonds_to_add )
    {
        frame.add_bond(bond.first[0], bond.first[1], bond.second);
    }
}

// TODO need to use std::string because read_multi_line can't return a string_view.
void mmCIFFormat::read_inline_property(const std::vector<string_view>& line_split, std::string& data)
{
    if ( line_split.size() == 2 )
    {
        data = line_split[1].to_string();
    }
    else
    {
        auto line = file_.readline();
        reader_meta_data_.line = line;

        if ( line[0] == ';' )
            data = read_multi_line();
        else
            data = line.to_string();
    }
}

void mmCIFFormat::fill_assembly()
{
    if ( category_map_.find("_pdbx_struct_assembly_gen") == category_map_.end() )
    { 
        return;
    }

    auto pdbx_struct_assembly_gen = category_map_["_pdbx_struct_assembly_gen"];
    file_.seekpos(pdbx_struct_assembly_gen.position_data_start);

    if ( pdbx_struct_assembly_gen.is_loop )
    {
       auto line = file_.readline();
       do
       {
           auto line_split = split(line, ' ');

           string_view assembly_id;
           string_view operation_expression;
           string_view targets_str;

           if ( line_split.size() == pdbx_struct_assembly_gen.property_map.size() )    // Single line data
           {
               assembly_id = line_split[pdbx_struct_assembly_gen.property_map["assembly_id"]];
               operation_expression = line_split[pdbx_struct_assembly_gen.property_map["oper_expression"]];
               targets_str = line_split[pdbx_struct_assembly_gen.property_map["asym_id_list"]];
           }
           else //Multiple line data
           {
               reader_meta_data_.line = line;

               size_t start_row_line = 0;
               size_t end_row_line = start_row_line + line_split.size() - 1;

               // Because we remove multiline symbols, we can't have a string_view.
               // In case of multiline, we save it in this variable and set a string_view on it.
               std::string multi_line;

               do
               {
                   end_row_line = start_row_line + line_split.size() - 1;
                   
                   if ( line[0] == ';' )
                   { 
                       multi_line = read_multi_line();
                       line_split[0] = multi_line;
                   }

                   auto row_index = pdbx_struct_assembly_gen.property_map["assembly_id"];
                   if ( start_row_line <= row_index && row_index <= end_row_line )
                   {
                       assembly_id = line_split[row_index - start_row_line];
                   }

                   row_index = pdbx_struct_assembly_gen.property_map["oper_expression"];
                   if ( start_row_line <= row_index && row_index <= end_row_line )
                   {
                       operation_expression = line_split[row_index - start_row_line];
                   }

                   row_index = pdbx_struct_assembly_gen.property_map["asym_id_list"];
                   if ( start_row_line <= row_index && row_index <= end_row_line )
                   {
                       targets_str = line_split[row_index - start_row_line];
                   }

                   start_row_line += line_split.size();

                   if ( start_row_line >= pdbx_struct_assembly_gen.property_map.size() )
                   {
                       break;
                   }

                   line = file_.readline();
                   reader_meta_data_.line = line;
                   line_split = split(line, ' ');
               } while ( file_.tellpos() < pdbx_struct_assembly_gen.position_data_end );
           }

           // Clean possible quote around oper_expression
           if ( operation_expression.size() > 2 && operation_expression[0] == '\'' && operation_expression[operation_expression.size() - 1] == '\'' )
           {
               operation_expression = operation_expression.substr(1, operation_expression.size() - 2);
           }

           std::set<std::string> targets = std::set<std::string>();
           fill_assembly_targets_vector(targets_str, targets);

           build_assembly_generators(assembly_id.to_string(), operation_expression, targets);
           line = file_.readline();
        } while ( file_.tellpos() < pdbx_struct_assembly_gen.position_data_end );
    }
    else
    {
        std::string assemblyGeneratorID;
        std::string operation_expression;
        std::set<std::string> targets;

        auto line = file_.readline();
        do
        {
            if ( line.empty() || line[0] == '#' )
                break;

            auto line_split = split(line, ' ');
            
            if (line_split[0] == "_pdbx_struct_assembly_gen.assembly_id")
            {
                read_inline_property(line_split, assemblyGeneratorID);
            }
            else if ( line_split[0] == "_pdbx_struct_assembly_gen.oper_expression" )
            {
                read_inline_property(line_split, operation_expression);
            }
            else if ( line_split[0] == "_pdbx_struct_assembly_gen.asym_id_list" )
            {
                std::string target_list;
                read_inline_property(line_split, target_list);
                fill_assembly_targets_vector(target_list, targets);
            }

            line = file_.readline();
        } while ( file_.tellpos() < pdbx_struct_assembly_gen.position_data_end );

        build_assembly_generators(assemblyGeneratorID, operation_expression, targets);
    }

    fill_assembly_operations();
}
void mmCIFFormat::fill_assembly_operations()
{
    auto struct_oper_list = category_map_["_pdbx_struct_oper_list"];
    file_.seekpos(struct_oper_list.position_data_start);

    auto struct_oper_list_map = struct_oper_list.property_map;

    auto operation_id = struct_oper_list_map.at("id");
    // auto operation_type = struct_oper_list_map.at("type");
    // auto operation_name = struct_oper_list_map.at("name");
    // auto symmetry_operation = struct_oper_list_map.at("symmetry_operation");

    auto matrix_1_1 = struct_oper_list_map.at("matrix[1][1]");
    auto matrix_1_2 = struct_oper_list_map.at("matrix[1][2]");
    auto matrix_1_3 = struct_oper_list_map.at("matrix[1][3]");
    auto matrix_2_1 = struct_oper_list_map.at("matrix[2][1]");
    auto matrix_2_2 = struct_oper_list_map.at("matrix[2][2]");
    auto matrix_2_3 = struct_oper_list_map.at("matrix[2][3]");
    auto matrix_3_1 = struct_oper_list_map.at("matrix[3][1]");
    auto matrix_3_2 = struct_oper_list_map.at("matrix[3][2]");
    auto matrix_3_3 = struct_oper_list_map.at("matrix[3][3]");

    auto vector_x = struct_oper_list_map.at("vector[1]");
    auto vector_y = struct_oper_list_map.at("vector[2]");
    auto vector_z = struct_oper_list_map.at("vector[3]");

    if ( struct_oper_list.is_loop )
    {
        std::vector<std::string> line_split = std::vector<std::string>();
        line_split.reserve(struct_oper_list_map.size());

        while (file_.tellpos() < struct_oper_list.position_data_end )
        {
            read_item_properties(struct_oper_list, line_split);

            if ( line_split.size() < struct_oper_list_map.size() )
            {
                throw format_error("_pdbx_struct_oper_list does not contains right number of properties ({} instead of {}) at position {} in '{}' ", 
                                   line_split.size(), struct_oper_list_map.size(), file_.tellpos(),
                                   file_.path());
            }

            auto rotation = Matrix3D(
                cif_to_double(line_split[matrix_1_1]), cif_to_double(line_split[matrix_1_2]), cif_to_double(line_split[matrix_1_3]),
                cif_to_double(line_split[matrix_2_1]), cif_to_double(line_split[matrix_2_2]), cif_to_double(line_split[matrix_2_3]),
                cif_to_double(line_split[matrix_3_1]), cif_to_double(line_split[matrix_3_2]), cif_to_double(line_split[matrix_3_3]));

            auto translation = Vector3D(cif_to_double(line_split[vector_x]),
                                        cif_to_double(line_split[vector_y]),
                                        cif_to_double(line_split[vector_z]));

            assembly_.assembly_operations.emplace(line_split[operation_id], AssemblyOperation(translation, rotation));
        }
    }
    else
    {
        string_view id_value;

        double matrix_1_1_value = 0.;
        double matrix_1_2_value = 0.;
        double matrix_1_3_value = 0.;
        double matrix_2_1_value = 0.;
        double matrix_2_2_value = 0.;
        double matrix_2_3_value = 0.;
        double matrix_3_1_value = 0.;
        double matrix_3_2_value = 0.;
        double matrix_3_3_value = 0.;

        double vector_x_value = 0.;
        double vector_y_value = 0.;
        double vector_z_value = 0.;

        while ( file_.tellpos() < struct_oper_list.position_data_end )
        {
            auto line = file_.readline();
            auto line_split = split(line, ' ');

            if ( line_split[0].ends_with(".id") )
            {
                id_value = line_split[1];
            }
            else if ( line_split[0].ends_with(".matrix[1][1]") )
            {
                matrix_1_1_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[1][2]") )
            {
                matrix_1_2_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[1][3]") )
            {
                matrix_1_3_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[2][1]") )
            {
                matrix_2_1_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[2][2]") )
            {
                matrix_2_2_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[2][3]") )
            {
                matrix_2_3_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[3][1]") )
            {
                matrix_3_1_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[3][2]") )
            {
                matrix_3_2_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".matrix[3][3]") )
            {
                matrix_3_3_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".vector[1]") )
            {
                vector_x_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".vector[2]") )
            {
                vector_y_value = cif_to_double(line_split[1].to_string());
            }
            else if ( line_split[0].ends_with(".vector[3]") )
            {
                vector_z_value = cif_to_double(line_split[1].to_string());
            }
        }

        auto rotation = Matrix3D(
            matrix_1_1_value, matrix_1_2_value, matrix_1_3_value,
            matrix_2_1_value, matrix_2_2_value, matrix_2_3_value, 
            matrix_3_1_value, matrix_3_2_value, matrix_3_3_value);

        auto translation = Vector3D(vector_x_value, vector_y_value, vector_z_value);

        assembly_.assembly_operations.emplace(id_value, AssemblyOperation(translation, rotation));
    }
}

void mmCIFFormat::write(const Frame& frame) {
    if (models_ == 0) {
        file_.print("# generated by Chemfiles\n");
        file_.print("#\n");
        auto lengths = frame.cell().lengths();
        file_.print("_cell.length_a {:#}\n", lengths[0]);
        file_.print("_cell.length_b {:#}\n", lengths[1]);
        file_.print("_cell.length_c {:#}\n", lengths[2]);
        auto angles = frame.cell().angles();
        file_.print("_cell.length_alpha {:#}\n", angles[0]);
        file_.print("_cell.length_beta  {:#}\n", angles[1]);
        file_.print("_cell.length_gamma {:#}\n", angles[2]);
        file_.print("#\n");
        file_.print("loop_\n");
        file_.print("_atom_site.group_PDB\n");
        file_.print("_atom_site.id\n");
        file_.print("_atom_site.type_symbol\n");
        file_.print("_atom_site.label_atom_id\n");
        file_.print("_atom_site.label_alt_id\n");
        file_.print("_atom_site.label_comp_id\n");
        file_.print("_atom_site.label_asym_id\n");
        file_.print("_atom_site.label_seq_id\n");
        file_.print("_atom_site.Cartn_x\n");
        file_.print("_atom_site.Cartn_y\n");
        file_.print("_atom_site.Cartn_z\n");
        file_.print("_atom_site.pdbx_formal_charge\n");
        file_.print("_atom_site.auth_asym_id\n");
        file_.print("_atom_site.auth_seq_id\n");
        file_.print("_atom_site.pdbx_PDB_model_num\n");
    }

    models_++;

    const auto& topology = frame.topology();
    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        ++atoms_;

        std::string compid = ".";
        std::string asymid = ".";
        std::string seq_id = ".";
        std::string auth_asymid = ".";
        std::string pdbgroup = "HETATM";

        const auto& residue = topology.residue_for_atom(i);
        if (residue) {
            compid = residue->name();

            if (residue->id()) {
                seq_id = std::to_string(*residue->id());
            } else {
                seq_id = "?";
            }

            asymid = residue->get<Property::STRING>("chainid").value_or("?");
            auth_asymid = residue->get<Property::STRING>("chainname").value_or(".");
            if (residue->get<Property::BOOL>("is_standard_pdb").value_or(false)) {
                pdbgroup = "ATOM  ";
            }
        }

        const auto& atom = frame[i];

        file_.print("{} {: <5} {: <2} {: <4} {} {: >3} {} {: >4} {:8.3f} {:8.3f} {:8.3f} {:#} {} {: >4} {}\n",
                pdbgroup, atoms_, atom.type(), atom.name(), ".", compid,
                asymid, seq_id, positions[i][0], positions[i][1], positions[i][2],
                atom.charge(), auth_asymid, seq_id, models_
        );
    }

}

std::string mmCIFFormat::read_multi_line()
{
    std::string result = "";
    auto line = reader_meta_data_.line;
    
    while ( !file_.eof() )
    {
        result += line.substr(1).to_string();

        line = file_.readline();
        trim(line);

        if ( line.empty() || line[0] == '#' || line == ";" )
        {
            break;
        }
    }

    return result;
}

void mmCIFFormat::fill_assembly_targets_vector(const string_view& targets_str, std::set<std::string>& targets)
{
    string_view target_group_sv;

    if ( targets_str[0] == '\'' && targets_str[targets_str.size() - 1] == '\'' )
    {
        target_group_sv = targets_str.substr(1, targets_str.size() - 2);
    }
    else
    {
        target_group_sv = targets_str;
    }

    auto target_vector = split(target_group_sv, ',');

    for ( auto target : target_vector )
        targets.emplace(target);
}
void mmCIFFormat::build_assembly_generators(const std::string& assembly_id,
                                            const string_view& operation_expression,
                                            const std::set<std::string>& targets)
{
    size_t special_symbol_position;
    if ( (special_symbol_position = operation_expression.find_last_of('(')) != std::string::npos )
    {
        auto closeParenthesisIndex = operation_expression.find_first_of(')', special_symbol_position);
        auto innerParenthesisExpressionSize = closeParenthesisIndex - special_symbol_position - 1;

        auto innerParenthesisExpression = operation_expression.substr(special_symbol_position+1, innerParenthesisExpressionSize);
        build_assembly_generators(assembly_id, innerParenthesisExpression, targets);

        auto lhs_expression = special_symbol_position == 0 ? "" : operation_expression.substr(0, special_symbol_position);
        auto rhs_expression = (closeParenthesisIndex == operation_expression.size() - 1) ? "" : operation_expression.substr(closeParenthesisIndex + 1);

        auto outer_expression = lhs_expression.to_string() + rhs_expression.to_string();

        if ( outer_expression.size() > 0 )
            build_assembly_generators(assembly_id, outer_expression, targets);
    }
    else if ( (special_symbol_position = operation_expression.find('-')) != std::string::npos )
    {
        auto firstTarget_str = operation_expression.substr(0, special_symbol_position);
        auto secondTarget_str = operation_expression.substr(special_symbol_position + 1);

        std::string firstTarget_string = firstTarget_str.to_string();
        std::string secondTarget_string = secondTarget_str.to_string();

        int firstTargetID = std::stoi(firstTarget_string);
        int secondTargetID = std::stoi(secondTarget_string);

        for ( int i = firstTargetID; i <= secondTargetID; ++i )
        {
            auto target_str = std::to_string(i);
            build_assembly_generators(assembly_id, target_str, targets);
        }
    }
    else if ( operation_expression.find(',') != std::string::npos )
    {
        auto operations = split(operation_expression, ',');

        for ( auto operation : operations )
        {
            build_assembly_generators(assembly_id, operation, targets);
        }
    }
    else
    {
        AssemblyGenerator assemblyGenerator = AssemblyGenerator();
        assemblyGenerator.assembly_id = assembly_id;
        assemblyGenerator.operations.emplace_back(operation_expression.to_string());
        for ( auto target : targets )
            assemblyGenerator.targets.emplace(target);
        assembly_.assembly_generators.emplace_back(assemblyGenerator);
    }
}

double cif_to_double(std::string line)
{
    line.erase(std::remove(line.begin(), line.end(), '('), line.end());
    line.erase(std::remove(line.begin(), line.end(), ')'), line.end());
    return parse<double>(line);
}

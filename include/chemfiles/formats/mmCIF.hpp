// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MMCIF_HPP
#define CHEMFILES_FORMAT_MMCIF_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>

#include "chemfiles/types.hpp"

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

#include "chemfiles/Residue.hpp"
#include "chemfiles/UnitCell.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// mmCIF Crystallographic Information Framework for MacroMolecules
/// reader and writer.
class mmCIFFormat final: public Format {

private:
    class ReaderMetaData
    {
    public:
        string_view line;
        size_t position;
    };

    class mmCIFCategoryHeader
    {
    public:
        std::map<std::string, size_t> property_map;
        bool is_loop;
        size_t position_data_start;
        size_t position_data_end;
    };

    class AssemblyOperation
    {
    public:
        AssemblyOperation() : AssemblyOperation(Vector3D(0, 0, 0), Matrix3D(1, 0, 0, 0, 1, 0 ,0, 0, 1))
        {}

        AssemblyOperation(const Vector3D& _translation, const Matrix3D& _rotation) : 
            translation(_translation), rotation(_rotation)
        {};

        Vector3D translation;
        Matrix3D rotation;
    };

    class AssemblyGenerator
    {
    public:
        std::string assembly_id;
        std::set<std::string> targets;
        std::vector<std::string> operations;
    };

    class Assembly
    {
    public:
        std::vector<AssemblyGenerator> assembly_generators;
        std::map<std::string, AssemblyOperation> assembly_operations;
    };

public:
    mmCIFFormat(std::string path, File::Mode mode, File::Compression compression) :
        file_(std::move(path), mode, compression), models_(0), atoms_(0) {
        init_();
    }

    mmCIFFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
        file_(std::move(memory), mode, compression), models_(0), atoms_(0) {
        init_();
    }

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;
    size_t nsteps() override;

private:
    /// Initialize important variables
    void init_();
    /// Underlying file representation
    TextFile file_;

    ReaderMetaData reader_meta_data_;
    std::map<std::string, mmCIFCategoryHeader> category_map_;

    Assembly assembly_;

    void init_atom_site();
    void addCategory(const std::string& category_name, bool is_loop);

    void read_atom_site(Frame& frame);
    void read_struct_oper_list(Frame& frame);

    void fill_assembly();
    void fill_assembly_operations();

    void apply_symmetry(Frame& frame, const std::string & assembly_id);

    void fill_loop_properties(const std::string& category_name, std::map<std::string, size_t>& property_map);

    void fill_targets_vector(const string_view& targets_str, std::set<string_view>& targets);
    void build_assembly_generators(const std::string& assembly_id, const string_view& operation_expression, const std::set<string_view>& targets);

    /// Map residues position, indexed by residue id and chainid. We use an indirection to keep the residue order (and don't sort them with the map id) 
    std::map<std::pair<std::string, int64_t>, size_t> map_residues_indexes;
    /// Store residues in this vector.
    std::vector<Residue> residues_;
    
    /// Storing the positions of all the steps in the file, so that we can
    /// just `seekpos` them instead of reading the whole step.
    std::vector<uint64_t> steps_positions_;
    /// The cell for all frames
    UnitCell cell_;
    /// Number of models written to the file.
    size_t models_;
    /// Number of atoms written to the file.
    size_t atoms_;
    /// Frame properties need to be stored
    std::string name_;
    /// The PDB icode, if any
    std::string pdb_idcode_;

    void read_item(const mmCIFCategoryHeader& header, std::vector<string_view>& property);
    void read_line_and_split(std::vector<string_view> & split) ;

    std::string read_multi_line();
};

template<> const FormatMetadata& format_metadata<mmCIFFormat>();

} // namespace chemfiles

#endif

// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_GRO_HPP
#define CHEMFILES_FORMAT_GRO_HPP

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/external/optional.hpp"
#include <cstdint>
#include <memory>
#include <string>

namespace chemfiles
{
	class Frame;
	class MemoryBuffer;
	class FormatMetadata;

	/// GRO file format reader and writer.
	class GROFormat final : public TextFormat
	{
	public:
		static const int64_t GRO_INDEX_MAX;

		GROFormat(std::string path, File::Mode mode, File::Compression compression) :
			TextFormat(std::move(path), mode, compression)
		{
		}

		GROFormat(std::shared_ptr<MemoryBuffer> memory, File::Mode mode, File::Compression compression) :
			TextFormat(std::move(memory), mode, compression)
		{
		}

		void			   read_next(Frame& frame) override;
		void			   write_next(const Frame& frame) override;
		optional<uint64_t> forward() override;

	private:
		static std::string get_chain_name_from_index(const int chain_index);
	};

	template<>
	const FormatMetadata& format_metadata<GROFormat>();

} // namespace chemfiles

#endif

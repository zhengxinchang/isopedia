1. Changed the profile format by merging the old and new profile files into a single format. Profile files from older versions (1.7.3 and earlier) are no longer compatible.
2. Added profile_info.json and profile_params.json alongside the profile file. They contain the parameters used to generate the profile file and a summary of it.
3. Replaced the 64-bit hash with a 128-bit hash to avoid hash collisions.
4. Removed the hash signature from PTIR to reduce its size.


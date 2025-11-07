-- =================================================================
--  Database Seed File
--  Generated from input_files document.
-- =================================================================

-- Drop existing objects --

DROP TABLE IF EXISTS trawl_report_catch CASCADE;
DROP TABLE IF EXISTS trawl_report_length CASCADE;
DROP TABLE IF EXISTS trawl_report_specimen CASCADE;
DROP TABLE IF EXISTS trawl_report_haul CASCADE;
DROP TABLE IF EXISTS species CASCADE;
DROP TYPE IF EXISTS sex_enum;

-- Create Custom Types --

CREATE TYPE sex_enum AS ENUM (
    'male',
    'female',
    'unsexed'
);

--  Create Lookup Tables --

CREATE TABLE trawl_report_haul (
    haul_num INTEGER PRIMARY KEY,
    cruise_name VARCHAR(100),
    haul_timestamp TIMESTAMPTZ DEFAULT NOW()
);

CREATE TABLE species (
    species_id INTEGER PRIMARY KEY,
    common_name VARCHAR(100),
    scientific_name VARCHAR(100)
);

-- Create Main Data Tables --

CREATE TABLE trawl_report_specimen (
    specimen_id SERIAL PRIMARY KEY,
    haul_num INTEGER NOT NULL REFERENCES trawl_report_haul(haul_num),
    species_id INTEGER NOT NULL REFERENCES species(species_id),

    sex sex_enum NOT NULL DEFAULT 'unsexed', 

    -- cm
    length DECIMAL(10, 2) CHECK (length > 0),
    
    -- kg
    weight DECIMAL(10, 3) CHECK (weight > 0),
    
    -- years
    age DECIMAL(5, 1) CHECK (age >= 0) 
);

CREATE TABLE trawl_report_length (
    length_dist_id SERIAL PRIMARY KEY,
    haul_num INTEGER NOT NULL REFERENCES trawl_report_haul(haul_num),
    species_id INTEGER NOT NULL REFERENCES species(species_id),

    sex sex_enum NOT NULL DEFAULT 'unsexed',

    -- cm
    length DECIMAL(10, 2) NOT NULL CHECK (length > 0),

    length_count INTEGER NOT NULL DEFAULT 0 CHECK (length_count >= 0)
);

CREATE TABLE trawl_report_catch (
    catch_id SERIAL PRIMARY KEY,
    haul_num INTEGER NOT NULL REFERENCES trawl_report_haul(haul_num),
    species_id INTEGER NOT NULL REFERENCES species(species_id),

    -- kg
    haul_weight DECIMAL(10, 3) NOT NULL CHECK (haul_weight >= 0),

    -- Ensure only one weight entry per haul/species
    UNIQUE(haul_num, species_id)
);

-- Insert Data --

INSERT INTO trawl_report_haul (haul_num, cruise_name) VALUES
(1, 'Summer 2024 Survey'),
(2, 'Summer 2024 Survey');

INSERT INTO species (species_id, common_name, scientific_name) VALUES
(306, 'Pacific Hake', 'Merluccius productus'),
(206, 'Walleye Pollock', 'Gadus chalcogrammus');

INSERT INTO trawl_report_specimen (haul_num, species_id, sex, length, weight, age) VALUES
(1, 306, 'male', 30.5, 0.450, 4.0),
(1, 306, 'male', 31.0, 0.465, 4.0),
(1, 306, 'female', 32.0, 0.510, 5.0),
(1, 306, 'unsexed', 15.2, NULL, 1.0),      -- NULL weight
(1, 206, 'female', 25.0, 0.300, 3.0),
(1, 206, 'female', 26.5, 0.320, 3.0),
(2, 306, 'male', 40.0, 0.600, 6.0),
(2, 306, 'female', 42.5, 0.650, 7.0),
(2, 306, 'unsexed', NULL, NULL, NULL); -- All info missing

INSERT INTO trawl_report_length (haul_num, species_id, sex, length, length_count) VALUES
(1, 306, 'male', 30.0, 2),
(1, 306, 'female', 32.0, 1),
(1, 306, 'unsexed', 15.0, 1),
(1, 206, 'female', 25.0, 2),
(2, 306, 'male', 40.0, 1),
(2, 306, 'female', 42.0, 1);

INSERT INTO trawl_report_catch (haul_num, species_id, haul_weight) VALUES
(1, 306, 120.500),
(1, 206, 75.200),
(2, 306, 250.000);
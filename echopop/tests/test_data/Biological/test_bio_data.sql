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
    ship INTEGER NOT NULL,
    survey INTEGER NOT NULL,
    haul_num INTEGER NOT NULL,
    haul_timestamp TIMESTAMPTZ DEFAULT NOW(),
    PRIMARY KEY (ship, survey, haul_num)
);

CREATE TABLE species (
    species_code INTEGER PRIMARY KEY,
    common_name VARCHAR(100),
    scientific_name VARCHAR(100)
);

-- Create Main Data Tables --

CREATE TABLE trawl_report_specimen (
    specimen_id SERIAL PRIMARY KEY,
    ship INTEGER NOT NULL,
    survey INTEGER NOT NULL,
    haul_num INTEGER NOT NULL,
    species_code INTEGER NOT NULL REFERENCES species(species_code),

    sex sex_enum NOT NULL DEFAULT 'unsexed',

    -- cm
    length DECIMAL(10, 2) CHECK (length > 0),

    -- kg
    weight DECIMAL(10, 3) CHECK (weight > 0),

    -- years
    age DECIMAL(5, 1) CHECK (age >= 0),
    FOREIGN KEY (ship, survey, haul_num) REFERENCES trawl_report_haul(ship, survey, haul_num)
);

CREATE TABLE trawl_report_length (
    length_dist_id SERIAL PRIMARY KEY,
    ship INTEGER NOT NULL,
    survey INTEGER NOT NULL,
    haul_num INTEGER NOT NULL,
    species_code INTEGER NOT NULL REFERENCES species(species_code),

    sex sex_enum NOT NULL DEFAULT 'unsexed',

    -- cm
    length DECIMAL(10, 2) NOT NULL CHECK (length > 0),

    length_count INTEGER NOT NULL DEFAULT 0 CHECK (length_count >= 0),
    FOREIGN KEY (ship, survey, haul_num) REFERENCES trawl_report_haul(ship, survey, haul_num)
);

CREATE TABLE trawl_report_catch (
    catch_id SERIAL PRIMARY KEY,
    ship INTEGER NOT NULL,
    survey INTEGER NOT NULL,
    haul_num INTEGER NOT NULL,
    species_code INTEGER NOT NULL REFERENCES species(species_code),

    -- kg
    haul_weight DECIMAL(10, 3) NOT NULL CHECK (haul_weight >= 0),

    -- Ensure only one weight entry per haul/species
    UNIQUE(ship, survey, haul_num, species_code),
    FOREIGN KEY (ship, survey, haul_num) REFERENCES trawl_report_haul(ship, survey, haul_num)
);

-- Insert Data --

INSERT INTO trawl_report_haul (ship, survey, haul_num) VALUES
(101, 2024, 1),
(101, 2024, 2),
(102, 2024, 1),
(101, 2025, 1);

INSERT INTO species (species_code, common_name, scientific_name) VALUES
(22500, 'Pacific Hake', 'Merluccius productus'),
(206, 'Walleye Pollock', 'Gadus chalcogrammus'),
(150, 'Dover Sole', 'Microstomus pacificus');

INSERT INTO trawl_report_specimen (ship, survey, haul_num, species_code, sex, length, weight, age) VALUES
(101, 2024, 1, 22500, 'male', 30.5, 0.450, 4.0),
(101, 2024, 1, 22500, 'male', 31.0, 0.465, 4.0),
(101, 2024, 1, 22500, 'female', 32.0, 0.510, 5.0),
(101, 2024, 1, 22500, 'unsexed', 15.2, NULL, 1.0),      -- NULL weight
(101, 2024, 1, 206, 'female', 25.0, 0.300, 3.0),
(101, 2024, 1, 206, 'female', 26.5, 0.320, 3.0),
(101, 2024, 2, 22500, 'male', 40.0, 0.600, 6.0),
(101, 2024, 2, 22500, 'female', 42.5, 0.650, 7.0),
(101, 2024, 2, 22500, 'unsexed', NULL, NULL, NULL), -- All info missing
(102, 2024, 1, 150, 'female', 45.0, 1.200, 10.0),
(102, 2024, 1, 150, 'male', 40.0, 0.950, 8.0),
(101, 2025, 1, 206, 'male', 35.0, 0.500, 5.0);

INSERT INTO trawl_report_length (ship, survey, haul_num, species_code, sex, length, length_count) VALUES
(101, 2024, 1, 22500, 'male', 30.0, 2),
(101, 2024, 1, 22500, 'female', 32.0, 1),
(101, 2024, 1, 22500, 'unsexed', 15.0, 1),
(101, 2024, 1, 206, 'female', 25.0, 2),
(101, 2024, 2, 22500, 'male', 40.0, 1),
(101, 2024, 2, 22500, 'female', 42.0, 1),
(102, 2024, 1, 150, 'unsexed', 42.0, 5),
(101, 2025, 1, 206, 'unsexed', 35.0, 3);

INSERT INTO trawl_report_catch (ship, survey, haul_num, species_code, haul_weight) VALUES
(101, 2024, 1, 22500, 120.500),
(101, 2024, 1, 206, 75.200),
(101, 2024, 2, 22500, 250.000),
(102, 2024, 1, 150, 50.000),
(101, 2025, 1, 206, 90.000);

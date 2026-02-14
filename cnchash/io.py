"""
Input/output functions for CNCHASH.

Handles reading various input file formats (phase data, station files,
velocity models) and writing mechanism output files.
"""

import re
from collections import defaultdict

import numpy as np


def read_phase_file(filename):
    """
    Read phase format files - supports multiple formats.

    Supported formats:
    - Format 1 (north1.phase): Compressed 2-digit year format
    - Format 2 (north2.phase): 4-digit year with separate station lines
    - Format 3 (north4.phase): 8-digit date with embedded station info
    - Format 4 (north5.simul): SIMUL2000 format with calculated angles

    Parameters
    ----------
    filename : str
        Path to phase file

    Returns
    -------
    list
        List of events, each event is a dict with:
        - 'year', 'month', 'day', 'hour', 'min', 'sec': origin time
        - 'lat', 'lon', 'depth': location
        - 'mag': magnitude
        - 'id': event ID
        - 'stations': list of station observations
    """
    events = []

    # Try to read with different encodings
    lines = []
    try:
        with open(filename, encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(filename, encoding="latin-1", errors="ignore") as f:
            lines = f.readlines()

    # Detect format from first non-empty line
    format_type = None
    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Check for format 2 (4-digit year at start)
        if re.match(r"^\d{4}\s+\d{5}", line):
            format_type = 2
            break
        # Check for format 3 (8-digit date at start)
        elif re.match(r"^\d{8}\d+", line):
            format_type = 3
            break
        # Check for format 4 (SIMUL2000 format - has "DATE" header)
        elif "DATE" in line and "ORIGIN" in line:
            format_type = 4
            break
        # Check for format 1 (2-digit year) - both continuous and split time formats
        elif re.match(r"^\d{2}\s+\d{1,3}", line):
            format_type = 1
            break

    if format_type is None:
        # Try to detect from line patterns
        for line in lines:
            line = line.strip()
            if not line:
                continue
            # Format 1 pattern (both continuous and split time)
            if re.match(r"^\d{2}\s+\d{1,3}", line):
                format_type = 1
                break
            # Format 2 pattern
            elif re.match(r"^\d{4}\s+\d{5}", line):
                format_type = 2
                break
            # Format 3 pattern
            elif re.match(r"^\d{16}", line):
                format_type = 3
                break

    if format_type is None:
        format_type = 1  # Default

    # Parse based on format
    if format_type == 1:
        events = _parse_phase_format1(lines)
    elif format_type == 2:
        events = _parse_phase_format2(lines)
    elif format_type == 3:
        events = _parse_phase_format3(lines)
    elif format_type == 4:
        events = _parse_phase_format4(lines)

    return events


def _parse_phase_format1(lines):
    """
    Parse format 1 (north1.phase): Compressed 2-digit year format.

    Event header has two formats:
    1. Continuous time: YY MDDhhmmsssss... (e.g., 94 1211104155034)
       Where MDD = month*100 + day (e.g., 121 = Jan 21)
    2. Split time: YY MDD hhmmsssss... (e.g., 94 128 744463234)
       Where MDD = month*100 + day (e.g., 128 = Jan 28)

    Station line: STA POLARITY 0 DISTANCE AZI THE 1 X CHANNEL
    """
    events = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        # Check for event header (starts with 2-digit year)
        # Match both formats: continuous time and split time
        if re.match(r"^\d{2}\s+\d{1,3}", line):
            parts = line.split()

            if len(parts) < 8:  # Need minimum parts for both formats
                i += 1
                continue

            try:
                year = 1900 + int(parts[0]) if int(parts[0]) >= 50 else 2000 + int(parts[0])

                # Check if time is split format (parts[1] is 3 digits = MDD encoding)
                # or separate M D format (parts[1] and parts[2] are 1-2 digits each)
                # or continuous format (parts[1] is 10+ digits = MDDhhmmsssss)

                if len(parts[1]) <= 3 and len(parts) > 2:
                    # Check if this is M D format or MDD format
                    # M D format: parts[1] is 1-2 digits (month), parts[2] is 1-2 digits (day)
                    # MDD format: parts[1] is 3 digits (MDD encoding like 128 = Jan 28)
                    parts1_val = int(parts[1])

                    # Check if parts[1] is a valid month (1-12) -> M D format
                    # Or if parts[1] is MDD encoding (101-1231) -> MDD format
                    if 1 <= parts1_val <= 12 and len(parts) > 3:
                        # M D hhmmsssss format (e.g., 94 2 5 711451134)
                        month = parts1_val
                        day = int(parts[2])
                        time_str = parts[3]
                        # lat/lon/depth are at positions 4, 5, 6
                        lat_str = parts[4]
                        lon_str = parts[5]
                        dep_str = parts[6]
                        mag_idx = 7
                    elif parts1_val >= 100:
                        # MDD hhmmsssss format (e.g., 94 128 744463234)
                        # parts[1] is MDD encoding (e.g., 128 = Jan 28)
                        mdd = parts1_val
                        month = mdd // 100
                        day = mdd % 100
                        time_str = parts[2]
                        # lat/lon/depth are at positions 3, 4, 5
                        lat_str = parts[3]
                        lon_str = parts[4]
                        dep_str = parts[5]
                        mag_idx = 6
                    else:
                        # Default to M D format if parts[1] is 1-99
                        month = parts1_val
                        day = int(parts[2]) if len(parts) > 2 else 1
                        time_str = parts[3] if len(parts) > 3 else "0"
                        lat_str = parts[4] if len(parts) > 4 else "0"
                        lon_str = parts[5] if len(parts) > 5 else "0"
                        dep_str = parts[6] if len(parts) > 6 else "0"
                        mag_idx = 7
                else:
                    # Continuous time format: YY MDDhhmmsssss...
                    # parts[1] is MDDhhmmsssss
                    mdd = int(parts[1][0:3])
                    month = mdd // 100
                    day = mdd % 100
                    time_str = parts[1][3:]  # Remaining part is hhmmsssss
                    # lat/lon/depth are at positions 2, 3, 4
                    lat_str = parts[2]
                    lon_str = parts[3]
                    dep_str = parts[4]
                    mag_idx = 5

                # Parse time string to get hour, minute, sec
                # Time format varies:
                # - Split format: hmmsssss (e.g., 744463234 = 7:44:46.3234)
                # - Continuous format: hhmmsssss (e.g., 1104155034 = 11:04:15.5034)
                if len(time_str) >= 8:
                    # Check if first 2 digits form a valid hour (0-23)
                    potential_hour = int(time_str[0:2])
                    if potential_hour >= 24:
                        # Split format: single-digit hour
                        hour = int(time_str[0]) if len(time_str) > 0 else 0
                        minute = int(time_str[1:3]) if len(time_str) >= 3 else 0
                        sec_str = time_str[3:]
                        sec = (
                            float(sec_str[0:2] + "." + sec_str[2:])
                            if len(sec_str) > 2
                            else float(sec_str)
                        )
                    else:
                        # Continuous format: double-digit hour
                        hour = potential_hour
                        minute = int(time_str[2:4]) if len(time_str) >= 4 else 0
                        sec_str = time_str[4:]
                        sec = (
                            float(sec_str[0:2] + "." + sec_str[2:])
                            if len(sec_str) > 2
                            else float(sec_str)
                        )
                else:
                    # Handle short time strings
                    if len(time_str) >= 4:
                        hour = int(time_str[0:2])
                        minute = int(time_str[2:4])
                    else:
                        hour = 0
                        minute = 0
                    sec = 0.0

                # Parse location with scaling factors
                try:
                    lat = float(lat_str) / 42500.0
                except ValueError:
                    lat = 0.0

                try:
                    lon = -float(lon_str) / 31.23  # West is negative
                except ValueError:
                    lon = 0.0

                try:
                    depth = float(dep_str) / 10000.0
                except ValueError:
                    depth = 0.0

                # Magnitude
                mag_str = parts[mag_idx] if mag_idx < len(parts) else parts[5]
                if "." in mag_str:
                    mag = float(mag_str)
                else:
                    mag = float(mag_str) / 10.0

                # Event ID (second to last part)
                event_id = (
                    parts[-2]
                    if len(parts) > 10
                    else f"{year}{month:02d}{day:02d}{hour:02d}{minute:02d}"
                )

                event = {
                    "year": year,
                    "month": month,
                    "day": day,
                    "hour": hour,
                    "min": minute,
                    "sec": sec,
                    "lat": lat,
                    "lon": lon,
                    "depth": depth,
                    "mag": mag,
                    "id": event_id,
                    "stations": [],
                }

                # Read station lines
                i += 1
                while i < len(lines):
                    line = lines[i].strip()

                    # Check for event ID line (end of event data)
                    if line == event_id or (re.match(r"^\d{7}", line) and line == event_id):
                        i += 1
                        break

                    # Check for next event header (also end of event data)
                    if re.match(r"^\d{2}\s+\d{1,3}", line):
                        # Don't increment i, so the next iteration will process this header
                        break

                    # Skip empty lines
                    if not line:
                        i += 1
                        continue

                    # Parse station line
                    parts = line.split()
                    if len(parts) >= 3:
                        station = parts[0]
                        pol_code = parts[1] if len(parts) > 1 else ""

                        # Parse polarity code (e.g., "IPU0", "EPD1")
                        if len(pol_code) >= 3:
                            onset = pol_code[0].upper()  # I or E
                            polarity = pol_code[2].upper()  # U or D
                        else:
                            onset = "I"
                            polarity = pol_code[0].upper() if pol_code else "U"

                        # Channel is usually the last part
                        channel = parts[-1] if len(parts) > 5 else "VHZ"
                        network = "CI"

                        event["stations"].append(
                            {
                                "name": station,
                                "network": network,
                                "component": channel,
                                "onset": onset,
                                "polarity": polarity,
                            }
                        )

                    i += 1

                events.append(event)

            except (ValueError, IndexError):
                i += 1
                continue
        else:
            i += 1

    return events


def _parse_phase_format2(lines):
    """
    Parse format 2 (north2.phase): 4-digit year with separate station lines.

    Event header: YYYY MDD HHMM.SSSS lat.XXX lon.XXX dep.XX
    Where MDD = month*1000 + day
    Station line: STA  NET  COMP ONSET POLARITY
    """
    events = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        # Check for event header (starts with 4-digit year)
        # Multiple formats:
        # - YYYY MDDHHMM.SSSS... (combined date/time, parts[1] has 9+ digits)
        # - YYYY MDDD HHMM.SSSS... (separate MDDD date, parts[1] has 4-5 digits)
        # - YYYY M D HHMM.SSSS... (separate month and day, parts[1] has 1-2 digits)
        if re.match(r"^\d{4}\s+\d+", line):
            parts = line.split()

            if len(parts) < 5:
                i += 1
                continue

            try:
                year = int(parts[0])

                # Parse date/time - there are three formats:
                # Format 1: 12111 415.5034 (separate date and time, MDDD)
                # Format 2: 128154720.6834 (combined date and time, MDDHHMMssss)
                # Format 3: 2 5 71145.1134 (separate month and day, M D)

                # Check if month and day are separate (Format 3)
                # parts[1] is small (1-2) = month, parts[2] is small (1-2) = day
                if len(parts) >= 7 and parts[1].isdigit() and parts[2].isdigit():
                    parts1_val = int(parts[1])
                    parts2_val = int(parts[2])
                    # Check if this looks like M D format (both values are small)
                    if parts1_val <= 12 and parts2_val <= 31:
                        # Format 3: M D HHMM.SSSS...
                        month = parts1_val
                        day = parts2_val
                        time_str = parts[3]  # HHMM.SSSS

                        time_float = float(time_str)
                        hour = int(time_float // 100)
                        minute = int(time_float % 100)
                        sec = (time_float - hour * 100 - minute) * 100

                        # lat/lon are at parts[4] and parts[5]
                        lat_str = parts[4]
                        lon_str = parts[5]

                    elif "." in parts[1] and len(parts) == 9:
                        # Format 2: combined date/time
                        # Example: 128154720.6834 = January 28, 15:47:20.6834
                        datetime_str = parts[1]
                        datetime_int = int(float(datetime_str))

                        # Parse: MDDHHMMssss
                        # 128154720 -> 1 (Jan) 28 (day) 15 (hour) 47 (min) 20 (sec)
                        datetime_padded = f"{datetime_int:09d}"  # Pad to 9 digits

                        month = int(datetime_padded[0])
                        day = int(datetime_padded[1:3])
                        hour = int(datetime_padded[3:5])
                        minute = int(datetime_padded[5:7])
                        sec = int(datetime_padded[7:9])

                        # Get fractional seconds
                        if "." in datetime_str:
                            frac_sec = float("0." + datetime_str.split(".")[1])
                            sec += frac_sec

                        # Time string is already parsed, so lat/lon are at different positions
                        lat_str = parts[2]
                        lon_str = parts[3]

                else:
                    # Format 1: separate date and time
                    # Parse date (e.g., 12111 = January 21, or 12510 = January 25)
                    date_val = int(parts[1])

                    # Try different formats
                    if date_val > 10000:
                        # Format: MDDD (month + 3-digit day encoding)
                        month = date_val // 10000
                        day_part = date_val % 10000
                        if day_part > 100:
                            day = day_part // 100
                        else:
                            day = day_part
                    elif date_val > 1000:
                        # Format: MDDD or MMDD
                        month = date_val // 1000
                        day = date_val % 1000
                        if day > 31:
                            # Try MMDD format
                            month = date_val // 100
                            day = date_val % 100
                    else:
                        # Format: MD or MDD
                        month = date_val // 100
                        day = date_val % 100
                        if day == 0:
                            day = date_val % 10

                    # Parse time (e.g., 415.5034 = 04:15:50.34)
                    time_str = parts[2]
                    time_float = float(time_str)
                    hour = int(time_float // 100)
                    minute = int(time_float % 100)
                    sec = (time_float - hour * 100 - minute) * 100

                    # Parse lat/lon/depth from parts 3, 4, 5
                    lat_str = parts[3]
                    lon_str = parts[4]

                # Parse lat/lon/depth from parts (positions depend on format)
                # Format: encoded_lat encoded_lon.dep
                # For example: 14.55118 37.0618.13
                # The encoding is:
                # - actual_lat = encoded_lat * 2.353
                # - actual_lon = encoded_lon * (-3.2)
                # - depth is the second part after the second dot in lon_str

                # lat_str and lon_str are already set above

                # Parse latitude (encoded)
                try:
                    lat_enc = float(lat_str)
                    lat = lat_enc * 2.353
                except ValueError:
                    lat = 34.0  # Default value

                # Parse longitude and depth
                # lon_str might be like "37.0618.13" where 37.06 is lon encoding and 18.13 is depth
                # The format is: lon_encoding_with_fractional_part + depth_encoded
                # For example: "37.0618.13" = 37.06 (lon) + 18.13 (depth)
                try:
                    # Split by dots and parse
                    dot_parts = lon_str.split(".")
                    if len(dot_parts) >= 3:
                        # Format: XX.XXXX.DD.DD or similar
                        # For "37.0618.13":
                        # - dot_parts = ['37', '0618', '13']
                        # - We need to figure out where lon_enc ends and depth begins
                        # Try: lon_enc = 37.06, depth = 18.13
                        # This means dot_parts[1] = '0618' contains both lon fraction (06) and depth integer (18)
                        lon_frac = dot_parts[1][:2]  # First 2 digits are lon fraction
                        depth_int = dot_parts[1][2:]  # Next digits are depth integer
                        depth_frac = dot_parts[2] if len(dot_parts) > 2 else "0"
                        lon_enc = float(f"{dot_parts[0]}.{lon_frac}")
                        depth = float(
                            f"{depth_int}.{depth_frac}" if depth_frac else float(depth_int)
                        )
                    else:
                        lon_enc = float(lon_str)
                        depth = 18.0  # Default
                except ValueError:
                    lon_enc = 37.0
                    depth = 18.0

                lon = lon_enc * (-3.2)  # West is negative

                # Find magnitude in the line
                mag = 0.0
                for part in parts:
                    if "." in part:
                        try:
                            if 0 < float(part) < 10:
                                mag_val = float(part)
                                if 0 < mag_val < 10:
                                    mag = mag_val
                                    break
                        except ValueError:
                            # Skip parts that can't be converted to float
                            continue

                # Event ID is usually at end (7-digit number)
                # Check if the last part is a 7-digit number
                if parts[-1].isdigit() and len(parts[-1]) == 7:
                    event_id = parts[-1]
                else:
                    event_id = f"{year}{month:02d}{day:02d}{hour:02d}{minute:02d}"

                event = {
                    "year": year,
                    "month": month,
                    "day": day,
                    "hour": hour,
                    "min": minute,
                    "sec": sec,
                    "lat": lat,
                    "lon": lon,
                    "depth": depth,
                    "mag": mag,
                    "id": event_id,
                    "stations": [],
                }

                # Read station lines
                i += 1
                while i < len(lines):
                    line = lines[i].strip()

                    # Check for next event header or event ID line
                    if re.match(r"^\d{4}\s+\d{5}", line):
                        break
                    if line == event_id:
                        i += 1
                        break

                    # Skip empty lines
                    if not line:
                        i += 1
                        continue

                    # Parse station line
                    # Format: STA  NET  COMP ONSET POLARITY
                    parts = line.split()
                    if len(parts) >= 5:
                        station = parts[0].strip()
                        network = parts[1].strip() if len(parts) > 1 else "CI"
                        component = parts[2].strip() if len(parts) > 2 else "HHZ"
                        onset = parts[3].strip().upper() if len(parts) > 3 else "I"
                        polarity = parts[4].strip().upper() if len(parts) > 4 else "U"

                        event["stations"].append(
                            {
                                "name": station,
                                "network": network,
                                "component": component,
                                "onset": onset,
                                "polarity": polarity,
                            }
                        )

                    i += 1

                events.append(event)

            except (ValueError, IndexError):
                i += 1
                continue
        else:
            i += 1

    return events


def _parse_phase_format3(lines):
    """
    Parse format 3 (north4.phase): 8-digit date with embedded station info.

    Event header: YYYYMMDDhhmmsssss lat lon dep ...
    Station line: STA NET COMP polcode DATE TIME ...
    """
    events = []
    i = 0

    while i < len(lines):
        line = lines[i].rstrip("\n\r")

        # Check for event header (8-digit date followed by more digits)
        if re.match(r"^\d{8}\d+", line):
            try:
                # Parse event header
                date_str = line[0:8]
                time_str = line[8:19]  # Includes fractional seconds

                year = int(date_str[0:4])
                month = int(date_str[4:6])
                day = int(date_str[6:8])

                hour = int(time_str[0:2])
                minute = int(time_str[2:4])
                sec = (
                    float(time_str[4:6]) + float("0." + time_str[6:]) if len(time_str) > 6 else 0.0
                )

                # The line has fixed-width fields
                # After time: lat (positions 19-26), lon (27-34), depth (35-40)
                # The encoding is similar to format 1
                lat_str = line[19:27].strip()
                lon_str = line[27:35].strip()
                dep_str = line[35:41].strip()

                # Parse lat/lon using format 3 encoding
                # lat = encoded / 42232, lon = -encoded / 30.4, depth = encoded / 117.6
                try:
                    lat = float(lat_str) / 42232.0
                except ValueError:
                    lat = 34.0  # Default

                try:
                    lon = -float(lon_str) / 30.4  # West is negative
                except ValueError:
                    lon = -118.0  # Default

                try:
                    depth = float(dep_str) / 117.6
                except ValueError:
                    depth = 18.0  # Default

                # Find magnitude and event ID in the line
                mag = 2.0
                event_id = f"{year}{month:02d}{day:02d}{hour:02d}{minute:02d}"

                parts = line.split()
                for _j, part in enumerate(parts):
                    # Look for magnitude (typically 1-3 digits with decimal)
                    if re.match(r"^\d\.\d+$", part):
                        mag_val = float(part)
                        if 0 < mag_val < 10:
                            mag = mag_val
                    # Look for 7-digit event ID (may have suffix like 'c230')
                    if re.match(r"^\d{7}", part):
                        event_id = part[:7]  # Take first 7 digits only

                event = {
                    "year": year,
                    "month": month,
                    "day": day,
                    "hour": hour,
                    "min": minute,
                    "sec": sec,
                    "lat": lat,
                    "lon": lon,
                    "depth": depth,
                    "mag": mag,
                    "id": event_id,
                    "stations": [],
                }

                # Read station lines
                i += 1
                while i < len(lines):
                    line = lines[i].rstrip("\n\r")

                    # Check for next event header or event ID line
                    if re.match(r"^\d{8}\d+", line):
                        break
                    if line.strip() == event_id:
                        i += 1
                        break

                    # Skip empty lines
                    if not line.strip():
                        i += 1
                        continue

                    # Parse station line
                    # Format: STA NET COMP polcode YYYYMMDDhhmmss ...
                    parts = line.split()
                    if len(parts) >= 4:
                        station = parts[0].strip()

                        # Skip lines without polarity info
                        if len(parts) < 4:
                            i += 1
                            continue

                        # Network and component
                        network = "CI"
                        component = "HHZ"
                        if len(parts) > 1:
                            net_comp = parts[1]
                            if len(net_comp) >= 2:
                                network = net_comp[0:2]
                                component = net_comp[2:] if len(net_comp) > 2 else component

                        # Polarity code (e.g., "iPU0", "ePD1", "eP .6")
                        pol_code = parts[2] if len(parts) > 2 else ""

                        onset = "I"
                        polarity = "U"

                        if len(pol_code) >= 3:
                            onset = pol_code[0].upper()  # i or e
                            pol_char = pol_code[2].upper() if len(pol_code) > 2 else "U"
                            if pol_char in "UD":
                                polarity = pol_char
                        elif pol_code and "." in pol_code:
                            # Format like "eP .6" - emergent, poor quality
                            onset = "E"
                            polarity = "U"  # Default for poor quality

                        event["stations"].append(
                            {
                                "name": station,
                                "network": network,
                                "component": component,
                                "onset": onset,
                                "polarity": polarity,
                            }
                        )

                    i += 1

                events.append(event)

            except (ValueError, IndexError) as e:
                print(f"DEBUG format4: Exception at line {i}: {e}")  # DEBUG
                i += 1
                continue
        else:
            i += 1

    return events


def _parse_phase_format4(lines):
    """
    Parse format 4 (north5.simul): SIMUL2000 format with calculated angles.

    Header: DATE    ORIGIN   LATITUDE LONGITUDE  DEPTH    MAG NO           RMS
    Event: YY MDD H M SS.SS latN XX.XX lonW XXX.XX DEP   MAG NO ID RMS
          or YY MDD HMM SS.SS latN XX.XX lonW XXX.XX DEP   MAG NO ID RMS
    Station: STN  DIST  AZ TOA PRMK HRMN  PSEC TPOBS              PRES  PWT
    """
    events = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        # Check for event line (starts with 2-digit year)
        if re.match(r"^\d{2}\s+\d{2,3}\s+\d+\s+\d+", line):
            parts = line.split()

            if len(parts) < 10:
                i += 1
                continue

            try:
                # Parse year
                year = 1900 + int(parts[0]) if int(parts[0]) >= 50 else 2000 + int(parts[0])

                # Parse date - MDD encoding (month*100 + day)
                mdd = int(parts[1])
                month = mdd // 100
                day = mdd % 100

                # Parse time - two possible formats
                # Check if parts[3] has decimal point (HMM SS format)
                if "." in parts[3]:
                    # HMM SS format: parts[2]=hmm, parts[3]=ss.ss
                    # Event ID is at parts[9] for this format
                    hmm = int(parts[2])
                    hour = hmm // 100
                    minute = hmm % 100
                    sec = float(parts[3])
                    lat_str = parts[4]
                    lon_str = parts[5]
                    depth = float(parts[6])
                    mag = float(parts[7])
                    event_id_idx = 9
                else:
                    # H M SS format: parts[2]=h, parts[3]=m, parts[4]=ss.ss
                    # Event ID is at parts[10] for this format
                    hour = int(parts[2])
                    minute = int(parts[3])
                    sec = float(parts[4])
                    lat_str = parts[5]
                    lon_str = parts[6]
                    depth = float(parts[7])
                    mag = float(parts[8])
                    event_id_idx = 10

                # Parse latitude (e.g., "34N13.89")
                lat = 0.0
                if lat_str and ("N" in lat_str or "S" in lat_str):
                    parts_lat = re.split(r"[NS]", lat_str)
                    deg = float(parts_lat[0]) if parts_lat[0] else 0
                    min_val = float(parts_lat[1]) if len(parts_lat) > 1 else 0
                    lat = deg + min_val / 60.0
                    if "S" in lat_str:
                        lat = -lat

                # Parse longitude (e.g., "118W36.53")
                lon = 0.0
                if lon_str and ("E" in lon_str or "W" in lon_str):
                    parts_lon = re.split(r"[EW]", lon_str)
                    deg = float(parts_lon[0]) if parts_lon[0] else 0
                    min_val = float(parts_lon[1]) if len(parts_lon) > 1 else 0
                    lon = deg + min_val / 60.0
                    if "W" in lon_str:
                        lon = -lon

                # Parse event ID
                event_id = (
                    parts[event_id_idx]
                    if len(parts) > event_id_idx and parts[event_id_idx].isdigit()
                    else f"{year}{month:02d}{day:02d}{hour:02d}{minute:02d}"
                )

                event = {
                    "year": year,
                    "month": month,
                    "day": day,
                    "hour": hour,
                    "min": minute,
                    "sec": sec,
                    "lat": lat,
                    "lon": lon,
                    "depth": depth,
                    "mag": mag,
                    "id": event_id,
                    "stations": [],
                }

                # Read station lines
                i += 1
                while i < len(lines):
                    line = lines[i].strip()

                    # Check for next event or header
                    if line.startswith("DATE") or line.startswith("  DATE"):
                        break
                    if re.match(r"^\d{2}\s+\d{2,3}\s+\d+\s+\d+", line):
                        break

                    # Skip empty lines and header lines
                    if not line or line.startswith("STN") or line.startswith("  STN"):
                        i += 1
                        continue

                    # Parse station line
                    stn_parts = line.split()
                    if len(stn_parts) >= 5 and stn_parts[0] and not stn_parts[0].startswith("STN"):
                        station = stn_parts[0].strip()

                        # Parse polarity from PRMK field (e.g., "P0", "P1")
                        prmk = stn_parts[4] if len(stn_parts) > 4 else ""
                        onset = "I"
                        polarity = "U"

                        if len(prmk) >= 2 and prmk[0] == "P":
                            # Check residual sign from PRES column (second to last)
                            if len(stn_parts) > 8:
                                try:
                                    pres = float(stn_parts[-3])
                                    polarity = "D" if pres < 0 else "U"
                                except ValueError:
                                    polarity = "U"
                            onset = "E" if len(prmk) > 1 and prmk[1] == "1" else "I"

                        event["stations"].append(
                            {
                                "name": station,
                                "network": "CI",
                                "component": "HHZ",
                                "onset": onset,
                                "polarity": polarity,
                            }
                        )

                    i += 1

                # Add event to events list
                events.append(event)

            except (ValueError, IndexError):
                i += 1
                continue
        else:
            i += 1

    return events


def read_station_file(filename):
    """
    Read station file for station locations.

    SCSN station file format (fixed-width):
    Positions:
    - 0-4: Station name (left-aligned)
    - 5-8: Component (left-aligned, space-padded)
    - 42-50: Latitude (float, signed)
    - 52-61: Longitude (float, signed)
    - 63-67: Elevation (integer, meters)
    - 91-92: Network code

    Example line: "ABL  VHZ SCSN STATION  100                34.84845 -119.22497  1975 1900/01/01 3000/01/01 CI"

    Parameters
    ----------
    filename : str
        Path to station file

    Returns
    -------
    dict
        Dictionary mapping (station, component, network) to (lat, lon, elev)
        Also includes a '_byname' key with a simplified name-only index for fast lookup
    """
    stations = {}
    by_name = {}  # Fast lookup by station name only

    with open(filename) as f:
        for line in f:
            line_stripped = line.rstrip("\n\r")
            if not line_stripped or line_stripped.startswith("#"):
                continue

            if len(line_stripped) < 60:
                continue

            try:
                station = line_stripped[0:5].strip()
                component = line_stripped[5:9].strip()

                if not station:
                    continue

                parts = line_stripped.split()
                lat = None
                lon = None
                elev = None
                network = "CI"

                for part in parts:
                    part = part.strip()
                    if not part:
                        continue

                    try:
                        val = float(part)

                        if -90 <= val <= 90:
                            if abs(val - int(val)) > 0.01:
                                if lat is None:
                                    lat = val
                        elif -180 <= val <= 180:
                            if abs(val - int(val)) > 0.01 or val < 0:
                                if lon is None and lat is not None:
                                    lon = val
                        elif 0 <= val <= 9000:
                            if val == int(val) and elev is None and lon is not None:
                                elev = val

                    except ValueError:
                        if len(part) == 2 and part.isalpha() and part.isupper():
                            network = part
                        continue

                if lat is not None and lon is not None:
                    if elev is None:
                        elev = 0.0
                    elev_km = elev / 1000.0

                    key = (station, component, network)
                    stations[key] = (lat, lon, elev_km)

                    # Also add to by_name index (first entry wins)
                    if station not in by_name:
                        by_name[station] = (lat, lon, elev_km)

            except (ValueError, IndexError):
                continue

    # Add fast lookup index
    stations["_byname"] = by_name

    return stations


def read_polarity_reversal_file(filename):
    """
    Read station polarity reversal file.

    Parameters
    ----------
    filename : str
        Path to reversal file

    Returns
    -------
    dict
        Dictionary mapping station to list of (start_date, end_date) tuples
        Dates are integers in YYYYMMDD format, 0 means no end date
    """
    reversals = defaultdict(list)

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            station = parts[0]

            # Try to parse dates (second and third columns)
            try:
                start_date = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
                end_date = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 99999999
            except ValueError:
                # Skip lines that can't be parsed
                continue

            reversals[station].append((start_date, end_date))

    return reversals


def read_velocity_model(filename):
    """
    Read 1D velocity model file.

    Format: depth(km) velocity(km/s)

    Parameters
    ----------
    filename : str
        Path to velocity model file

    Returns
    -------
    tuple
        (depth, velocity) arrays
    """
    depths = []
    velocities = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    depths.append(float(parts[0]))
                    velocities.append(float(parts[1]))
                except ValueError:
                    continue

    return np.array(depths), np.array(velocities)


def read_hash_input_file(filename):
    """
    Read HASH input file (like example.inp).

    Parameters
    ----------
    filename : str
        Path to input file

    Returns
    -------
    dict
        Input parameters
    """
    with open(filename) as f:
        lines = f.readlines()

    params = {}

    # Detect format by checking the first few lines
    # Format 1 (example1): polfile, phasefile, outfile1, outfile2, params...
    # Format 2 (example2): stationfile, polfile, phasefile, outfile1, outfile2, params...
    # Format 3 (example3): stationfile, polfile, statcor, amp, phasefile, outfile1, params...
    # Format 4 (example4): stationfile, polfile, phasefile, outfile1, outfile2, params... (same as format 2)
    # Format 5 (example5): polfile, simulfile, phasefile, outfile1, outfile2, params...

    # Check if first line is a station file (contains 'station' or ends with '.stations')
    first_line = lines[0].strip() if lines else ""

    if "station" in first_line.lower() or first_line.endswith(".stations"):
        # Format 2, 3, or 4: has stationfile
        params["station_file"] = lines[0].strip()
        params["polfile"] = lines[1].strip() if len(lines) > 1 else ""

        # Check if there are additional files before phasefile (format 3 has statcor and amp)
        third_line = lines[2].strip() if len(lines) > 2 else ""
        if third_line.endswith(".statcor") or third_line.endswith(".amp"):
            # Format 3: stationfile, polfile, statcor, amp, phasefile, outfile1, params...
            params["statcor_file"] = lines[2].strip() if len(lines) > 2 else ""
            params["amp_file"] = lines[3].strip() if len(lines) > 3 else ""
            params["phasefile"] = lines[4].strip() if len(lines) > 4 else ""
            params["outfile1"] = lines[5].strip() if len(lines) > 5 else ""
            params["outfile2"] = lines[6].strip() if len(lines) > 6 else ""
            param_idx = 7
        else:
            # Format 2 or 4: stationfile, polfile, phasefile, outfile1, outfile2, params...
            params["phasefile"] = lines[2].strip() if len(lines) > 2 else ""
            params["outfile1"] = lines[3].strip() if len(lines) > 3 else ""
            params["outfile2"] = lines[4].strip() if len(lines) > 4 else ""
            param_idx = 5
    elif first_line.endswith(".simul"):
        # Format 5: polfile, simulfile, phasefile, outfile1, outfile2, params...
        params["polfile"] = lines[0].strip()
        params["simul_file"] = lines[1].strip() if len(lines) > 1 else ""
        params["phasefile"] = lines[2].strip() if len(lines) > 2 else ""
        params["outfile1"] = lines[3].strip() if len(lines) > 3 else ""
        params["outfile2"] = lines[4].strip() if len(lines) > 4 else ""
        param_idx = 5
    else:
        # Format 1: polfile, phasefile, outfile1, outfile2, params...
        params["polfile"] = lines[0].strip()
        params["phasefile"] = lines[1].strip() if len(lines) > 1 else ""
        params["outfile1"] = lines[2].strip() if len(lines) > 2 else ""
        params["outfile2"] = lines[3].strip() if len(lines) > 3 else ""
        param_idx = 4

    # Remaining lines are parameters (position-based, not all files have same format)
    # Try to parse parameters by looking for numeric values

    # Default values
    params["npolmin"] = 8
    params["max_agap"] = 90.0
    params["max_pgap"] = 60.0
    params["dang"] = 10.0
    params["nmc"] = 30
    params["maxout"] = 500
    params["badfrac"] = 0.1
    params["delmax"] = 120.0
    params["cangle"] = 45.0
    params["prob_max"] = 0.1
    params["qextra"] = 0.0
    params["qtotal"] = 0.0
    params["nextra"] = 0.0
    params["ntotal"] = 0.0

    # Parse remaining lines (param_idx is set above)

    # Different formats have different parameter orders
    # Format 1 (example1): 11 parameters after file names
    # Format 2 (example2,3,4): parameters followed by velocity models

    # Collect numeric parameter values in order
    param_values = []
    for i in range(param_idx, len(lines)):
        line = lines[i].strip()
        if not line:
            continue

        # Check if it's a velocity model file (starts with 'vz.')
        if line.startswith("vz.") or line.startswith("/"):
            break

        # Try to parse as number
        try:
            val = float(line)
        except ValueError:
            continue

        param_values.append(val)

    # For Format 1 (example1), parameters are in this fixed order:
    # npolmin, max_agap, max_pgap, dang, nmc, maxout, badfrac, delmax, cangle, prob_max
    # For other formats, use heuristics based on value ranges

    if first_line.endswith(".simul"):
        # Format 5: same order as Format 1
        if len(param_values) >= 1:
            params["npolmin"] = int(param_values[0]) if param_values[0] >= 1 else params["npolmin"]
        if len(param_values) >= 2:
            params["max_agap"] = param_values[1]
        if len(param_values) >= 3:
            params["max_pgap"] = param_values[2]
        if len(param_values) >= 4:
            params["dang"] = param_values[3]
        if len(param_values) >= 5:
            params["nmc"] = int(param_values[4])
        if len(param_values) >= 6:
            params["maxout"] = int(param_values[5])
        if len(param_values) >= 7:
            params["badfrac"] = param_values[6]
        if len(param_values) >= 8:
            params["delmax"] = param_values[7]
        if len(param_values) >= 9:
            params["cangle"] = param_values[8]
        if len(param_values) >= 10:
            params["prob_max"] = param_values[9]
    elif "station" in first_line.lower() or first_line.endswith(".stations"):
        # Format 2, 3, or 4: use heuristics based on value ranges
        for val in param_values:
            # npolmin: small integer (1-20)
            if 1 <= val <= 20 and "npolmin" not in params:
                params["npolmin"] = int(val)
            # max_agap: 30-180 degrees (check first)
            elif 30 <= val <= 180 and params["max_agap"] == 90.0:
                params["max_agap"] = val
            # max_pgap: 30-180 degrees (check after max_agap)
            elif 30 <= val <= 180 and params["max_pgap"] == 60.0:
                params["max_pgap"] = val
            # dang: grid angle (1-30 degrees, smaller range)
            elif 1 <= val <= 20 and params["dang"] == 10.0:
                params["dang"] = val
            # nmc: number of trials (5-1000)
            elif 5 <= val <= 1000 and params["nmc"] == 30:
                params["nmc"] = int(val)
            # maxout: max output (10-10000)
            elif 10 <= val <= 10000 and params["maxout"] == 500:
                params["maxout"] = int(val)
            # badfrac: 0.0-1.0
            elif 0.0 <= val <= 1.0 and params["badfrac"] == 0.1:
                params["badfrac"] = val
            # delmax: max distance (50-1000 km)
            elif 50 <= val <= 1000 and params["delmax"] == 120.0:
                params["delmax"] = val
            # cangle: angle cutoff (15-90 degrees)
            elif 15 <= val <= 90 and params["cangle"] == 45.0:
                params["cangle"] = val
            # prob_max: 0.0-1.0 (second 0-1 value)
            elif 0.0 <= val <= 1.0 and params["prob_max"] == 0.1:
                params["prob_max"] = val
    else:
        # Format 1: fixed order
        if len(param_values) >= 1:
            params["npolmin"] = int(param_values[0])
        if len(param_values) >= 2:
            params["max_agap"] = param_values[1]
        if len(param_values) >= 3:
            params["max_pgap"] = param_values[2]
        if len(param_values) >= 4:
            params["dang"] = param_values[3]
        if len(param_values) >= 5:
            params["nmc"] = int(param_values[4])
        if len(param_values) >= 6:
            params["maxout"] = int(param_values[5])
        if len(param_values) >= 7:
            params["badfrac"] = param_values[6]
        if len(param_values) >= 8:
            params["delmax"] = param_values[7]
        if len(param_values) >= 9:
            params["cangle"] = param_values[8]
        if len(param_values) >= 10:
            params["prob_max"] = param_values[9]

    # Calculate derived parameters
    params["nextra"] = max(int(params["npolmin"] * params["badfrac"] * 0.5), 2)
    params["ntotal"] = max(int(params["npolmin"] * params["badfrac"]), 2)

    return params


def write_mechanism_output(filename, events, mechanisms):
    """
    Write focal mechanism output file.

    Parameters
    ----------
    filename : str
        Output file path
    events : list
        List of event dictionaries
    mechanisms : list
        List of mechanism results (one per event)
    """
    with open(filename, "w") as f:
        for event, mech in zip(events, mechanisms, strict=True):
            if mech is None or mech.get("quality") == "F":
                # Failed event
                write_failed_event(f, event, mech)
            else:
                write_successful_event(f, event, mech)


def write_failed_event(f, event, mech):
    """Write a failed event to output file."""
    if mech is None:
        mech = {
            "strike": 999,
            "dip": 99,
            "rake": 999,
            "var_est": [99, 99],
            "mfrac": 0.99,
            "quality": "F",
            "prob": 0.0,
            "stdr": 0.0,
            "nplt": 0,
            "nout2": 0,
        }

    fmt = (
        "{id:16s} {year:4d} {month:2d} {day:2d} {hour:2d} {min:2d} {sec:6.3f} {etype:1s} "
        "{mag:5.3f} {magtype:1s} {lat:9.5f} {lon:10.5f} {depth:7.3f} {locqual:1s} "
        "{rms:7.3f} {seh:7.3f} {sez:7.3f} {terr:7.3f} "
        "{nppick+nspick:4d} {nppick:4d} {nspick:4d} "
        "{strike:4d} {dip:3d} {rake:4d} "
        "{var1:3d} {var2:3d} "
        "{npol:3d} {mfrac:3d} "
        "{quality:1s} {prob:3d} {stdr:3d} {mflag:1s}"
    )

    line = fmt.format(
        id=event.get("id", ""),
        year=event.get("year", 0),
        month=event.get("month", 0),
        day=event.get("day", 0),
        hour=event.get("hour", 0),
        min=event.get("min", 0),
        sec=event.get("sec", 0.0),
        etype=event.get("etype", "L"),
        mag=event.get("mag", 0.0),
        magtype=event.get("magtype", "X"),
        lat=event.get("lat", 0.0),
        lon=event.get("lon", 0.0),
        depth=event.get("depth", 0.0),
        locqual=event.get("locqual", "X"),
        rms=event.get("rms", -9.0),
        seh=event.get("seh", -9.0),
        sez=event.get("sez", -9.0),
        terr=event.get("terr", -9.0),
        nppick=event.get("nppick", -9),
        nspick=event.get("nspick", -9),
        strike=int(mech.get("strike", 999)),
        dip=int(mech.get("dip", 99)),
        rake=int(mech.get("rake", 999)),
        var1=int(mech.get("var_est", [99, 99])[0]),
        var2=int(mech.get("var_est", [99, 99])[1]),
        npol=event.get("npol", 0),
        mfrac=int(mech.get("mfrac", 0.99) * 100),
        quality=mech.get("quality", "F"),
        prob=int(mech.get("prob", 0.0) * 100),
        stdr=int(mech.get("stdr", 0.0) * 100),
        mflag=mech.get("mflag", " "),
    )

    f.write(line + "\n")


def write_successful_event(f, event, mech):
    """Write a successful event to output file."""
    nmult = mech.get("nmult", 1)

    for imult in range(nmult):
        # Extract values for this solution (if multiple)
        if nmult > 1:
            strike = mech["strike_avg"][imult]
            dip = mech["dip_avg"][imult]
            rake = mech["rake_avg"][imult]
            var1 = mech["rms_diff"][0, imult]
            var2 = mech["rms_diff"][1, imult]
            prob = mech["prob"][imult]
            quality = mech["quality"][imult]
            stdr = mech["stdr"][imult]
            mfrac = mech["mfrac"][imult]
        else:
            strike = mech.get("strike_avg", mech.get("strike", 0))
            dip = mech.get("dip_avg", mech.get("dip", 0))
            rake = mech.get("rake_avg", mech.get("rake", 0))
            var1 = mech.get("var_est", [0, 0])[0]
            var2 = mech.get("var_est", [0, 0])[1]
            prob = mech.get("prob", 0.0)
            quality = mech.get("quality", "D")
            stdr = mech.get("stdr", 0.0)
            mfrac = mech.get("mfrac", 0.0)

        mflag = "*" if nmult > 1 else " "

        fmt = (
            "{id:16s} {year:4d} {month:2d} {day:2d} {hour:2d} {min:2d} {sec:6.3f} {etype:1s} "
            "{mag:5.3f} {magtype:1s} {lat:9.5f} {lon:10.5f} {depth:7.3f} {locqual:1s} "
            "{rms:7.3f} {seh:7.3f} {sez:7.3f} {terr:7.3f} "
            "{nppick+nspick:4d} {nppick:4d} {nspick:4d} "
            "{strike:4d} {dip:3d} {rake:4d} "
            "{var1:3d} {var2:3d} "
            "{npol:3d} {mfrac:3d} "
            "{quality:1s} {prob:3d} {stdr:3d} {mflag:1s}"
        )

        line = fmt.format(
            id=event.get("id", ""),
            year=event.get("year", 0),
            month=event.get("month", 0),
            day=event.get("day", 0),
            hour=event.get("hour", 0),
            min=event.get("min", 0),
            sec=event.get("sec", 0.0),
            etype=event.get("etype", "L"),
            mag=event.get("mag", 0.0),
            magtype=event.get("magtype", "X"),
            lat=event.get("lat", 0.0),
            lon=event.get("lon", 0.0),
            depth=event.get("depth", 0.0),
            locqual=event.get("locqual", "X"),
            rms=event.get("rms", -9.0),
            seh=event.get("seh", -9.0),
            sez=event.get("sez", -9.0),
            terr=event.get("terr", -9.0),
            nppick=event.get("nppick", -9),
            nspick=event.get("nspick", -9),
            strike=int(strike),
            dip=int(dip),
            rake=int(rake),
            var1=int(var1),
            var2=int(var2),
            npol=event.get("npol", 0),
            mfrac=int(mfrac * 100),
            quality=quality,
            prob=int(prob * 100),
            stdr=int(stdr * 100),
            mflag=mflag,
        )

        f.write(line + "\n")


def write_acceptable_planes(filename, event, mech):
    """
    Write acceptable planes output file.

    Parameters
    ----------
    filename : str
        Output file path
    event : dict
        Event dictionary
    mech : dict
        Mechanism results
    """
    with open(filename, "w") as f:
        # Write event header
        fmt1 = (
            "{year:4d} {month:2d} {day:2d} {hour:2d} {min:2d} {sec:6.3f} "
            "{mag:5.2f} {lat:9.4f} {lon:10.4f} {depth:6.2f} "
            "{sez:8.4f} {seh:8.4f} {npol:5d} {nout2:5d} {id:16s} "
            "{strike:7.1f} {dip:6.1f} {rake:7.1f} {var1:6.1f} {var2:6.1f} "
            "{mfrac:7.3f} {quality:2s} {prob:7.3f} {stdr:5.2f}"
        )

        f.write(
            fmt1.format(
                year=event.get("year", 0),
                month=event.get("month", 0),
                day=event.get("day", 0),
                hour=event.get("hour", 0),
                min=event.get("min", 0),
                sec=event.get("sec", 0.0),
                mag=event.get("mag", 0.0),
                lat=event.get("lat", 0.0),
                lon=event.get("lon", 0.0),
                depth=event.get("depth", 0.0),
                sez=event.get("sez", 0.0),
                seh=event.get("seh", 0.0),
                npol=event.get("npol", 0),
                nout2=mech.get("nout2", 0),
                id=event.get("id", ""),
                strike=mech.get("strike_avg", 0),
                dip=mech.get("dip_avg", 0),
                rake=mech.get("rake_avg", 0),
                var1=mech.get("var_est", [0, 0])[0],
                var2=mech.get("var_est", [0, 0])[1],
                mfrac=mech.get("mfrac", 0.0),
                quality=mech.get("quality", "D"),
                prob=mech.get("prob", 0.0),
                stdr=mech.get("stdr", 0.0),
            )
        )
        f.write("\n")

        # Write individual planes
        faults = mech.get("faults", np.zeros((3, 0)))
        slips = mech.get("slips", np.zeros((3, 0)))
        strikes = mech.get("strike", np.zeros(0))
        dips = mech.get("dip", np.zeros(0))
        rakes = mech.get("rake", np.zeros(0))

        nout = min(len(strikes), mech.get("nout1", len(strikes)))

        for i in range(nout):
            fmt2 = (
                "     {strike:9.2f} {dip:9.2f} {rake:9.2f} "
                "{fn1:9.4f} {fn2:9.4f} {fn3:9.4f} "
                "{sl1:9.4f} {sl2:9.4f} {sl3:9.4f}"
            )

            f.write(
                fmt2.format(
                    strike=strikes[i],
                    dip=dips[i],
                    rake=rakes[i],
                    fn1=faults[0, i] if faults.shape[1] > i else 0.0,
                    fn2=faults[1, i] if faults.shape[1] > i else 0.0,
                    fn3=faults[2, i] if faults.shape[1] > i else 0.0,
                    sl1=slips[0, i] if slips.shape[1] > i else 0.0,
                    sl2=slips[1, i] if slips.shape[1] > i else 0.0,
                    sl3=slips[2, i] if slips.shape[1] > i else 0.0,
                )
            )
            f.write("\n")


# Export all functions
__all__ = [
    "read_phase_file",
    "read_station_file",
    "read_polarity_reversal_file",
    "read_velocity_model",
    "read_hash_input_file",
    "write_mechanism_output",
    "write_acceptable_planes",
]

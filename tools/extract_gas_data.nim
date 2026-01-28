## Extract gas cross-section data from Magboltz Fortran source
## Parses DATA statements and metadata from GAS subroutines
## Outputs structured JSON files for each gas

import std/[strutils, strscans, tables, json, os, re, parseutils, sequtils, algorithm, sets]

type
  DataArray = object
    name: string
    values: seq[float]

  GasData = object
    id: int
    name: string
    version: string
    mass_amu: float
    polarizability: float
    am2: float  # Born-Bethe parameter
    c: float    # Born-Bethe parameter
    nIon: int
    nIn: int    # Number of inelastic (excitation) channels
    nAtt: int   # Number of attachment channels
    arrays: Table[string, seq[float]]
    ionThresholds: seq[float]  # EION values
    eoby: seq[float]           # Opal-Beaty parameters
    excitationThresholds: seq[float]
    excitationNames: seq[string]

const
  GasNames = [
    (1, "cf4"),
    (2, "argon"),
    (3, "helium4"),
    (4, "helium3"),
    (5, "neon"),
    (6, "krypton"),
    (7, "xenon"),
    (8, "methane"),
    (9, "ethane"),
    (10, "propane"),
    (11, "isobutane"),
    (12, "co2"),
    (13, "neopentane"),
    (14, "h2o"),
    (15, "oxygen"),
    (16, "nitrogen"),
    (17, "nitric_oxide"),
    (18, "nitrous_oxide"),
    (19, "ethene"),
    (20, "acetylene"),
    (21, "hydrogen"),
    (22, "deuterium"),
    (23, "carbon_monoxide"),
    (24, "methylal"),
    (25, "dme"),
    (26, "reid_step"),
    (27, "maxwell"),
    (28, "reid_ramp"),
    (29, "c2f6"),
    (30, "sf6"),
    (31, "nh3"),
    (32, "propene"),
    (33, "cyclopropane"),
    (34, "methanol"),
    (35, "ethanol"),
    (36, "isopropanol"),
    (37, "npropanol"),
    (38, "cesium"),
    (39, "fluorine"),
    (40, "cs2"),
    (41, "cosh2"),
    (42, "cosf2"),
    (43, "cos2f10"),
    (44, "tma"),
    (45, "tea"),
    (46, "para_hydrogen"),
    (47, "ortho_deuterium"),
    (48, "brf3"),
    (49, "c3f8"),
    (50, "ozone"),
    (51, "mercury"),
    (52, "h2s"),
    (53, "n_butane"),
    (54, "sih4"),
    (55, "geh4"),
    (56, "sif4"),
    (57, "nf3"),
    (58, "c2h2f4"),
    (59, "ccl4"),
  ]

proc parseFortranFloat(s: string): tuple[value: float, ok: bool] =
  ## Parse Fortran floating point with D exponent notation
  ## Returns (value, true) on success, (0.0, false) on failure
  var cleaned = s.strip()

  # Quick check - if it doesn't look like a number, fail fast
  if cleaned.len == 0:
    return (0.0, false)

  # Must start with digit, '.', '-', or '+'
  if cleaned[0] notin {'0'..'9', '.', '-', '+'}:
    return (0.0, false)

  # Handle Fortran D notation (e.g., 1.0D5, 1.0D-5)
  cleaned = cleaned.replace("D", "e").replace("d", "e")
  # Handle cases like .5 -> 0.5
  if cleaned.len > 0 and cleaned[0] == '.':
    cleaned = "0" & cleaned
  elif cleaned.len > 1 and cleaned[0] == '-' and cleaned[1] == '.':
    cleaned = "-0" & cleaned[1..^1]
  elif cleaned.len > 1 and cleaned[0] == '+' and cleaned[1] == '.':
    cleaned = "0" & cleaned[1..^1]

  try:
    result = (parseFloat(cleaned), true)
  except ValueError:
    result = (0.0, false)

proc extractDataValues(content: string): seq[float] =
  ## Extract numeric values from DATA statement content
  result = @[]
  var i = 0
  var current = ""
  var depth = 0

  while i < content.len:
    let c = content[i]
    if c == '(':
      depth.inc
      i.inc
      continue
    elif c == ')':
      depth.dec
      i.inc
      continue

    if depth == 0:
      if c in {',', '/', '\n', ' ', '\t'}:
        if current.strip().len > 0:
          # Check for repeat notation like "27*.0115" meaning 27 copies of .0115
          let trimmed = current.strip()
          if '*' in trimmed:
            let parts = trimmed.split('*')
            if parts.len == 2:
              try:
                let count = parseInt(parts[0].strip())
                let (value, ok) = parseFortranFloat(parts[1])
                if ok:
                  for j in 0..<count:
                    result.add(value)
              except: discard
            else:
              let (value, ok) = parseFortranFloat(trimmed)
              if ok:
                result.add(value)
          else:
            let (value, ok) = parseFortranFloat(trimmed)
            if ok:
              result.add(value)
          current = ""
        i.inc
        continue
      else:
        current.add(c)
    i.inc

  # Handle last value
  if current.strip().len > 0:
    let trimmed = current.strip()
    if '*' in trimmed:
      let parts = trimmed.split('*')
      if parts.len == 2:
        try:
          let count = parseInt(parts[0].strip())
          let (value, ok) = parseFortranFloat(parts[1])
          if ok:
            for j in 0..<count:
              result.add(value)
        except: discard
    else:
      let (value, ok) = parseFortranFloat(trimmed)
      if ok:
        result.add(value)

proc findSubroutineRange(lines: seq[string], gasNum: int): (int, int) =
  ## Find the line range for a GAS subroutine
  let pattern = "SUBROUTINE GAS" & $gasNum & "("
  var startLine = -1
  var endLine = -1

  for i, line in lines:
    if pattern in line.toUpperAscii():
      startLine = i
    elif startLine >= 0 and "SUBROUTINE GAS" in line.toUpperAscii() and "(" in line:
      # Found next subroutine
      endLine = i - 1
      break
    elif startLine >= 0 and "END" == line.strip().toUpperAscii():
      endLine = i
      break

  if startLine >= 0 and endLine < 0:
    endLine = lines.len - 1

  result = (startLine, endLine)

proc isComment(line: string): bool =
  ## Check if a Fortran line is a comment
  ## Comments have C, c, *, or ! in column 1
  if line.len == 0:
    return false
  return line[0] in {'C', 'c', '*', '!'}

proc extractDataStatements(lines: seq[string], startLine, endLine: int): Table[string, seq[float]] =
  ## Extract all DATA statements from a subroutine
  result = initTable[string, seq[float]]()
  var i = startLine

  while i <= endLine:
    let line = lines[i]

    # Skip comment lines
    if isComment(line):
      i.inc
      continue

    let trimmed = line.strip()
    let upper = trimmed.toUpperAscii()

    # Look for DATA statements (not in comments)
    if upper.startsWith("DATA "):
      # Extract array name
      let afterData = trimmed[5..^1].strip()
      let slashPos = afterData.find('/')
      if slashPos >= 0:
        let arrayName = afterData[0..<slashPos].strip().toUpperAscii()
        var content = afterData[slashPos + 1..^1]

        # Check if this DATA statement continues on following lines
        # In Fortran fixed format, continuation lines have:
        # - Columns 1-5: blank or label
        # - Column 6: continuation character (non-blank)
        # - Column 7+: continued statement
        # For DATA statements, continuation often looks like "     /values"
        var j = i + 1
        while j <= endLine:
          let nextLine = lines[j]

          # Skip comment lines
          if isComment(nextLine):
            j.inc
            continue

          # Check for continuation (non-blank in column 6 or line starting with whitespace then /)
          let trimmed = nextLine.strip()
          if nextLine.len >= 6 and nextLine[5] notin {' ', '0', '\t'}:
            # Standard Fortran continuation - column 6 has continuation char
            var contContent = nextLine
            if contContent.len > 6:
              contContent = contContent[6..^1]
            content.add(" " & contContent.strip())
            j.inc
          elif trimmed.startsWith("/") and trimmed.len > 1:
            # Alternative: line starts with / after whitespace
            content.add(" " & trimmed)
            j.inc
          else:
            break

        # Remove trailing slash
        if content.endsWith("/"):
          content = content[0..^2]

        # Also handle cases where there are multiple values sections
        let values = extractDataValues(content)
        if values.len > 0:
          result[arrayName] = values
    i.inc

proc extractMetadata(lines: seq[string], startLine, endLine: int): tuple[
  name: string,
  mass: float,
  apol: float,
  am2: float,
  c: float,
  nion: int,
  nin: int,
  natt: int,
  eion: seq[float],
  eoby: seq[float]
] =
  ## Extract metadata from assignment statements
  var
    name = ""
    mass = 0.0
    apol = 0.0
    am2 = 0.0
    c = 0.0
    nion = 0
    nin = 0
    natt = 0
    eion: seq[float] = @[]
    eoby: seq[float] = @[]

  for i in startLine..endLine:
    let line = lines[i]
    let trimmed = line.strip()
    let upper = trimmed.toUpperAscii()

    # NAME assignment
    if "NAME=" in upper or "NAME =" in upper:
      let quoteStart = trimmed.find('\'')
      if quoteStart >= 0:
        let quoteEnd = trimmed.find('\'', quoteStart + 1)
        if quoteEnd > quoteStart:
          name = trimmed[quoteStart + 1..<quoteEnd].strip()

    # E(2) gives mass info: E(2) = 2.0*EMASS/(mass*AMU)
    # So mass = 2.0*EMASS/E(2)/AMU, but we can extract from E(2)=2.0*EMASS/(X*AMU)
    if "E(2)=" in upper or "E(2) =" in upper:
      # Try to extract mass from pattern like "2.0*EMASS/(39.948*AMU)"
      let massMatch = trimmed.find(re"([0-9.]+)\s*\*\s*AMU")
      if massMatch >= 0:
        var m: array[1, string]
        if trimmed.match(re".*?([0-9.]+)\s*\*\s*AMU.*", m):
          try:
            mass = parseFloat(m[0])
          except: discard

    # APOL assignment
    if upper.startsWith("APOL=") or upper.startsWith("APOL ="):
      let eqPos = trimmed.find('=')
      if eqPos >= 0:
        let valStr = trimmed[eqPos + 1..^1].strip()
        let (val, ok) = parseFortranFloat(valStr)
        if ok:
          apol = val

    # AM2 assignment
    if upper.startsWith("AM2=") or upper.startsWith("AM2 ="):
      let eqPos = trimmed.find('=')
      if eqPos >= 0:
        let valStr = trimmed[eqPos + 1..^1].strip()
        let (val, ok) = parseFortranFloat(valStr)
        if ok:
          am2 = val

    # C assignment (careful - C is common)
    if (upper.startsWith("C=") or upper.startsWith("C =")) and not upper.startsWith("C ") and upper[0] == 'C':
      let eqPos = trimmed.find('=')
      if eqPos >= 0 and eqPos < 3:  # C= should be near start
        let valStr = trimmed[eqPos + 1..^1].strip()
        let (val, ok) = parseFortranFloat(valStr)
        if ok:
          c = val

    # NION assignment
    if upper.startsWith("NION=") or upper.startsWith("NION ="):
      let eqPos = trimmed.find('=')
      if eqPos >= 0:
        let valStr = trimmed[eqPos + 1..^1].strip()
        try:
          nion = parseInt(valStr)
        except: discard

    # NIN assignment
    if upper.startsWith("NIN=") or upper.startsWith("NIN ="):
      let eqPos = trimmed.find('=')
      if eqPos >= 0:
        let valStr = trimmed[eqPos + 1..^1].strip()
        try:
          nin = parseInt(valStr)
        except: discard

    # NATT assignment
    if upper.startsWith("NATT=") or upper.startsWith("NATT ="):
      let eqPos = trimmed.find('=')
      if eqPos >= 0:
        let valStr = trimmed[eqPos + 1..^1].strip()
        try:
          natt = parseInt(valStr)
        except: discard

    # EION(n) assignments
    if "EION(" in upper:
      let m = trimmed.find(re"EION\s*\(\s*(\d+)\s*\)\s*=\s*([0-9.DEde+-]+)")
      if m >= 0:
        var matches: array[2, string]
        if trimmed.toUpperAscii().match(re".*EION\s*\(\s*(\d+)\s*\)\s*=\s*([0-9.DEde+-]+).*", matches):
          try:
            let idx = parseInt(matches[0]) - 1
            let (val, ok) = parseFortranFloat(matches[1])
            if ok:
              while eion.len <= idx:
                eion.add(0.0)
              eion[idx] = val
          except: discard

    # EOBY(n) assignments
    if "EOBY(" in upper:
      var matches: array[2, string]
      if trimmed.toUpperAscii().match(re".*EOBY\s*\(\s*(\d+)\s*\)\s*=\s*([0-9.DEde+-]+).*", matches):
        try:
          let idx = parseInt(matches[0]) - 1
          let (val, ok) = parseFortranFloat(matches[1])
          if ok:
            while eoby.len <= idx:
              eoby.add(0.0)
            eoby[idx] = val
        except: discard

  result = (name, mass, apol, am2, c, nion, nin, natt, eion, eoby)

proc categorizeArrays(arrays: Table[string, seq[float]]): JsonNode =
  ## Categorize extracted arrays into elastic, ionization, excitation etc.
  result = newJObject()
  var usedArrays = initHashSet[string]()

  # Elastic scattering data
  var elastic = newJObject()

  # Momentum transfer - rare gas format
  if "XEN" in arrays and "YSEC" in arrays:
    elastic["momentum_transfer"] = %*{
      "energy_eV": arrays["XEN"],
      "cross_section_A2": arrays["YSEC"]
    }
    usedArrays.incl("XEN")
    usedArrays.incl("YSEC")

  # Momentum transfer - molecular format (YELM)
  if "XEN" in arrays and "YELM" in arrays:
    elastic["momentum_transfer"] = %*{
      "energy_eV": arrays["XEN"],
      "cross_section_A2": arrays["YELM"]
    }
    usedArrays.incl("XEN")
    usedArrays.incl("YELM")

  # Total elastic - rare gas format (YEL)
  if "XEN" in arrays and "YEL" in arrays:
    elastic["total"] = %*{
      "energy_eV": arrays["XEN"],
      "cross_section_A2": arrays["YEL"]
    }
    usedArrays.incl("YEL")

  # Total elastic - molecular format (YELT)
  if "XEN" in arrays and "YELT" in arrays:
    elastic["total"] = %*{
      "energy_eV": arrays["XEN"],
      "cross_section_A2": arrays["YELT"]
    }
    usedArrays.incl("YELT")

  # Angular distribution epsilon
  if "XEPS" in arrays and "YEPS" in arrays:
    elastic["angular_epsilon"] = %*{
      "energy_eV": arrays["XEPS"],
      "epsilon": arrays["YEPS"]
    }
    usedArrays.incl("XEPS")
    usedArrays.incl("YEPS")

  if elastic.len > 0:
    result["elastic"] = elastic

  # Ionization data
  var ionization = newJObject()
  var ionChannels = newJArray()

  # Ground state ionization (XENI/YENI, YENC, YEN1) - rare gas format
  if "XENI" in arrays:
    usedArrays.incl("XENI")
    if "YENI" in arrays:
      ionChannels.add(%*{
        "name": "gross",
        "energy_eV": arrays["XENI"],
        "cross_section_A2": arrays["YENI"]
      })
      usedArrays.incl("YENI")
    if "YENC" in arrays:
      ionChannels.add(%*{
        "name": "counting",
        "energy_eV": arrays["XENI"],
        "cross_section_A2": arrays["YENC"]
      })
      usedArrays.incl("YENC")
    if "YEN1" in arrays:
      ionChannels.add(%*{
        "name": "charge_state_1",
        "energy_eV": arrays["XENI"],
        "cross_section_A2": arrays["YEN1"]
      })
      usedArrays.incl("YEN1")

  # Molecular ionization format (XION/YION, YINC)
  if "XION" in arrays:
    usedArrays.incl("XION")
    if "YION" in arrays:
      ionChannels.add(%*{
        "name": "total",
        "energy_eV": arrays["XION"],
        "cross_section_A2": arrays["YION"]
      })
      usedArrays.incl("YION")
    if "YINC" in arrays:
      ionChannels.add(%*{
        "name": "counting",
        "energy_eV": arrays["XION"],
        "cross_section_A2": arrays["YINC"]
      })
      usedArrays.incl("YINC")

    # Dissociative ionization fragments (YIN1-YIN20)
    for i in 1..20:
      let yName = "YIN" & $i
      if yName in arrays:
        ionChannels.add(%*{
          "name": "fragment_" & $i,
          "energy_eV": arrays["XION"],
          "cross_section_A2": arrays[yName]
        })
        usedArrays.incl(yName)

  # Fragment ionization (XINF/YINF patterns for CH4 etc)
  if "XINF" in arrays and "YINF" in arrays:
    ionChannels.add(%*{
      "name": "fragment_total",
      "energy_eV": arrays["XINF"],
      "cross_section_A2": arrays["YINF"]
    })
    usedArrays.incl("XINF")
    usedArrays.incl("YINF")

  for i in 1..10:
    let xName = "XINF" & $i
    let yName = "YINF" & $i
    if xName in arrays and yName in arrays:
      ionChannels.add(%*{
        "name": "fragment_" & $i,
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # Pair production (XINPP/YINPP)
  if "XINPP" in arrays and "YINPP" in arrays:
    ionChannels.add(%*{
      "name": "pair_production",
      "energy_eV": arrays["XINPP"],
      "cross_section_A2": arrays["YINPP"]
    })
    usedArrays.incl("XINPP")
    usedArrays.incl("YINPP")

  # Higher charge states (XEN2/YEN2, XEN3/YEN3)
  for i in 2..10:
    let xName = "XEN" & $i
    let yName = "YEN" & $i
    if xName in arrays and yName in arrays:
      ionChannels.add(%*{
        "name": "charge_state_" & $i,
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # Shell ionizations
  var shells = newJArray()

  # K-shell (multiple formats)
  if "XKSH" in arrays and "YKSH" in arrays:
    shells.add(%*{"name": "K", "energy_eV": arrays["XKSH"], "cross_section_A2": arrays["YKSH"]})
    usedArrays.incl("XKSH")
    usedArrays.incl("YKSH")

  # K-shell for specific atoms (Carbon, Fluorine, Sulfur, etc)
  for suffix in ["C", "F", "S", "N", "O"]:
    let xName = "XKSH" & suffix
    let yName = "YKSH" & suffix
    if xName in arrays and yName in arrays:
      shells.add(%*{"name": "K_" & suffix, "energy_eV": arrays[xName], "cross_section_A2": arrays[yName]})
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # L-shells (multiple formats)
  for i in 1..3:
    for suffix in ["S", "SH"]:
      let xName = "XL" & $i & suffix
      let yName = "YL" & $i & suffix
      if xName in arrays and yName in arrays:
        shells.add(%*{"name": "L" & $i, "energy_eV": arrays[xName], "cross_section_A2": arrays[yName]})
        usedArrays.incl(xName)
        usedArrays.incl(yName)

  # M-shells
  for i in 1..5:
    for suffix in ["S", "SH"]:
      let xName = "XM" & $i & suffix
      let yName = "YM" & $i & suffix
      if xName in arrays and yName in arrays:
        shells.add(%*{"name": "M" & $i, "energy_eV": arrays[xName], "cross_section_A2": arrays[yName]})
        usedArrays.incl(xName)
        usedArrays.incl(yName)

  # N-shells (for heavy elements)
  for i in 1..7:
    let xName = "XN" & $i & "S"
    let yName = "YN" & $i & "S"
    if xName in arrays and yName in arrays:
      shells.add(%*{"name": "N" & $i, "energy_eV": arrays[xName], "cross_section_A2": arrays[yName]})
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  if ionChannels.len > 0:
    ionization["channels"] = ionChannels
  if shells.len > 0:
    ionization["shells"] = shells

  if ionization.len > 0:
    result["ionization"] = ionization

  # Excitation data - look for paired X/Y arrays
  var excitation = newJArray()

  # Rare gas excitation states
  let excitationPatterns = [
    "1S5", "1S4", "1S3", "1S2",
    "2S5", "2S4", "2S3", "2S2",
    "3S1PPPP", "3S1PPP", "3S1PP", "3S1P", "3S1",
    "2P10", "2P9", "2P8", "2P7", "2P6", "2P5", "2P4", "2P3", "2P2", "2P1",
    "3P10", "3P9", "3P8", "3P7", "3P6", "3P5", "3P4", "3P3", "3P2", "3P1",
    "3D6", "3D5", "3D4P", "3D4", "3D3", "3D2", "3D1PP", "3D1P", "3D1",
    "4D6", "4D5", "4D4P", "4D4", "4D3", "4D2", "4D1PP", "4D1P", "4D1",
    "4F", "5F",
    "SUM", "SUMP", "SUMD", "SUMF"
  ]

  for pattern in excitationPatterns:
    let xName = "X" & pattern
    let yName = "Y" & pattern
    if xName in arrays and yName in arrays:
      excitation.add(%*{
        "name": pattern,
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # Vibrational modes - multiple naming conventions
  # XVBV1-4/YVBV1-4 (bending), XVBH1-2/YVBH1-2 (stretching)
  for prefix in ["VBV", "VBH", "VIB"]:
    for i in 1..10:
      let xName = "X" & prefix & $i
      let yName = "Y" & prefix & $i
      if xName in arrays and yName in arrays:
        excitation.add(%*{
          "name": prefix & $i,
          "type": "vibrational",
          "energy_eV": arrays[xName],
          "cross_section_A2": arrays[yName]
        })
        usedArrays.incl(xName)
        usedArrays.incl(yName)

  # SF6 style vibrational (XV1V1, XV2V1, etc)
  for i in 1..10:
    let xName = "XV" & $i & "V1"
    let yName = "YV" & $i & "V1"
    if xName in arrays and yName in arrays:
      excitation.add(%*{
        "name": "V" & $i & "V1",
        "type": "vibrational",
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # Triplet states (XTR1-3/YTR1-3)
  for i in 1..10:
    let xName = "XTR" & $i
    let yName = "YTR" & $i
    if xName in arrays and yName in arrays:
      excitation.add(%*{
        "name": "triplet_" & $i,
        "type": "triplet",
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # Molecular dissociation products
  let dissocPatterns = [
    ("XDET", "YDET", "detachment"),
    ("XCHB", "YCHB", "CH_B_band"),
    ("XCHD", "YCHD", "CH_dissociation"),
    ("XHAL", "YHAL", "H_alpha"),
    ("XHBE", "YHBE", "H_beta")
  ]
  for (xName, yName, name) in dissocPatterns:
    if xName in arrays and yName in arrays:
      excitation.add(%*{
        "name": name,
        "type": "dissociation",
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  if excitation.len > 0:
    result["excitation"] = excitation

  # Attachment data
  var attachment = newJArray()
  if "XATT" in arrays and "YATT" in arrays:
    attachment.add(%*{
      "name": "total",
      "energy_eV": arrays["XATT"],
      "cross_section_A2": arrays["YATT"]
    })
    usedArrays.incl("XATT")
    usedArrays.incl("YATT")

  # Multiple attachment channels (XATT1/YATT1 format)
  for i in 1..8:
    let xName = "XATT" & $i
    let yName = "YATT" & $i
    if xName in arrays and yName in arrays:
      attachment.add(%*{
        "name": "channel_" & $i,
        "energy_eV": arrays[xName],
        "cross_section_A2": arrays[yName]
      })
      usedArrays.incl(xName)
      usedArrays.incl(yName)

  # SF6 style attachment (YAT1-7 with XION as energy scale)
  if "XION" in arrays:
    for i in 1..10:
      let yName = "YAT" & $i
      if yName in arrays:
        attachment.add(%*{
          "name": "channel_" & $i,
          "energy_eV": arrays["XION"],
          "cross_section_A2": arrays[yName]
        })
        usedArrays.incl(yName)

  if attachment.len > 0:
    result["attachment"] = attachment

  # Store any remaining unmatched arrays
  var otherArrays = newJObject()
  for name, values in arrays:
    if name notin usedArrays:
      otherArrays[name] = %values

  if otherArrays.len > 0:
    result["other_arrays"] = otherArrays

proc extractGas(lines: seq[string], gasNum: int): JsonNode =
  ## Extract all data for a single gas
  let (startLine, endLine) = findSubroutineRange(lines, gasNum)

  if startLine < 0:
    echo "  Warning: GAS", gasNum, " subroutine not found"
    return nil

  echo "  Found GAS", gasNum, " at lines ", startLine + 1, "-", endLine + 1

  # Extract DATA statements
  let arrays = extractDataStatements(lines, startLine, endLine)
  echo "    Extracted ", arrays.len, " arrays"

  # Extract metadata
  let (name, mass, apol, am2, c, nion, nin, natt, eion, eoby) =
    extractMetadata(lines, startLine, endLine)

  # Build JSON output
  result = newJObject()
  result["id"] = %gasNum
  result["name"] = %name
  result["mass_amu"] = %mass

  if apol > 0:
    result["polarizability"] = %apol

  if am2 > 0 or c > 0:
    result["born_bethe"] = %*{
      "AM2": am2,
      "C": c
    }

  result["channels"] = %*{
    "ionization": nion,
    "excitation": nin,
    "attachment": natt
  }

  if eion.len > 0:
    result["ionization_thresholds_eV"] = %eion

  if eoby.len > 0:
    result["opal_beaty"] = %eoby

  # Categorize and add cross-section data
  let categorized = categorizeArrays(arrays)
  for key, val in categorized:
    result[key] = val

proc main() =
  let magboltzPath = "/home/basti/src/Magboltz/magboltz-11.19.f"
  let outputDir = "/home/basti/CastData/ExternCode/NimBoltz/gas_data"

  echo "Reading Magboltz source: ", magboltzPath
  let content = readFile(magboltzPath)
  let lines = content.splitLines()
  echo "  Total lines: ", lines.len

  # Find all GAS subroutines
  var subroutineNums: seq[int] = @[]
  for i, line in lines:
    let upper = line.toUpperAscii()
    if "SUBROUTINE GAS" in upper:
      # Extract number from SUBROUTINE GASnn(
      var num: int
      let match = upper.find("SUBROUTINE GAS")
      if match >= 0:
        let afterGas = upper[match + 14..^1]
        var numStr = ""
        for c in afterGas:
          if c.isDigit:
            numStr.add(c)
          else:
            break
        if numStr.len > 0:
          try:
            num = parseInt(numStr)
            if num > 0 and num notin subroutineNums:
              subroutineNums.add(num)
          except: discard

  subroutineNums.sort()
  echo "\nFound ", subroutineNums.len, " GAS subroutines: ", subroutineNums

  # Create output directory
  createDir(outputDir)

  # Extract each gas
  for gasNum in subroutineNums:
    echo "\nProcessing GAS", gasNum, "..."

    let gasData = extractGas(lines, gasNum)
    if gasData != nil:
      # Find file name
      var fileName = "gas" & align($gasNum, 2, '0')
      for (id, name) in GasNames:
        if id == gasNum:
          fileName = fileName & "_" & name
          break
      fileName = fileName & ".json"

      let outputPath = outputDir / fileName
      writeFile(outputPath, gasData.pretty())
      echo "  Wrote: ", outputPath

when isMainModule:
  main()

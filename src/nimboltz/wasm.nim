## NimBoltz in a web page (Emscripten): a `Setup` as JSON in, the `Results`
## as JSON out. Built by `web/build.sh`; `web/worker.js` calls it in a Web
## Worker, so a run doesn't freeze the page, and passes Magboltz's printout
## (stdout) on line by line.
##
## Setup: {"mixture": [[gas number, %], ...], "temperature": K,
## "pressure": Torr, "efield": V/cm, "bfield": T, "angle": degrees,
## "collisions": n, "efinal": eV (0: found), "penning": bool,
## "thermal": bool, "seed": n, "report": bool (printout to stdout)}; all
## but "mixture" and "efield" optional, with the defaults of `Setup`.
##
## Results: {"gases", "efinal", "collisions", "nullCollisions",
## "velocity": {x, y, z}, "diffusion": {xx, yy, zz, xy, xz, yz},
## "longitudinal", "transverse", "meanEnergy", "ionisation",
## "attachment", "characteristicEnergy": {longitudinal, transverse},
## "spread": {...}, "frequencies", "processes", "spectrum": {energy,
## fraction}, "warnings"}, each value {"value", "error", "unit"}; or
## {"error": message} when the run stops.

import std/[json, options]
import api, gases, vec3
from drift import Tensor

proc measurement[T](m: Measurement[T]): JsonNode =
  %*{"value": float64(m.value), "error": float64(m.error), "unit": unitOf(m.value)}

proc quantity[T](x: T): JsonNode = %*{"value": float64(x), "unit": unitOf(x)}

proc vector[T](v: Vector[T]): JsonNode =
  %*{"x": measurement(v.x), "y": measurement(v.y), "z": measurement(v.z)}

proc tensor[T](d: Tensor[T]): JsonNode =
  %*{"xx": measurement(d.xx), "yy": measurement(d.yy), "zz": measurement(d.zz),
     "xy": measurement(d.xy), "xz": measurement(d.xz), "yz": measurement(d.yz)}

proc toJson(r: Results): JsonNode =
  result = %*{
    "gases": r.gases, "efield": quantity(r.efield), "efinal": quantity(r.efinal),
    "collisions": r.collisions, "nullCollisions": r.nullCollisions,
    "velocity": vector(r.velocity), "diffusion": tensor(r.diffusion),
    "longitudinal": measurement(r.longitudinal), "transverse": measurement(r.transverse),
    "meanEnergy": measurement(r.meanEnergy),
    "ionisation": measurement(r.ionisation), "attachment": measurement(r.attachment),
    "characteristicEnergy": {"longitudinal": measurement(r.characteristicEnergy(r.longitudinal)),
                             "transverse": measurement(r.characteristicEnergy(r.transverse))},
    "spread": {"longitudinal": measurement(r.spread(r.longitudinal)),
               "transverse": measurement(r.spread(r.transverse))},
    "frequencies": {"total": quantity(r.frequencies.total),
                    "elastic": quantity(r.frequencies.elastic),
                    "ionisation": quantity(r.frequencies.ionisation),
                    "attachment": quantity(r.frequencies.attachment),
                    "inelastic": quantity(r.frequencies.inelastic)},
    "warnings": r.warnings,
    "analyses": {"steadyState": r.steadyState.isSome,
                 "pulsedTownsend": r.pulsedTownsend.isSome,
                 "timeOfFlight": r.timeOfFlight.isSome}}
  var processes = newJArray()
  for p in r.processes:
    processes.add %*{"gas": p.gas, "description": p.description,
                     "frequency": measurement(p.frequency)}
  result["processes"] = processes
  var energy, fraction = newJArray()
  var last = 0
  for k, (e, f) in r.spectrum:
    if f > 0.0: last = k
  for k in 0 .. last:
    energy.add %float64(r.spectrum[k].energy)
    fraction.add %r.spectrum[k].fraction
  result["spectrum"] = %*{"energy": energy, "fraction": fraction, "unit": "eV"}

proc setup(j: JsonNode): Setup =
  result = Setup(efield: j["efield"].getFloat.V•cm⁻¹)
  for g in j["mixture"]:
    let id = g[0].getInt
    if id notin ord(Gas.low) .. ord(Gas.high):
      raise newException(MagboltzError, "no gas number " & $id)
    result.mixture.add (Gas(id), g[1].getFloat)
  if j.hasKey("temperature"): result.temperature = j["temperature"].getFloat.K
  if j.hasKey("pressure"): result.pressure = j["pressure"].getFloat.Torr
  if j.hasKey("bfield"): result.bfield = j["bfield"].getFloat.T
  if j.hasKey("angle"): result.angle = j["angle"].getFloat
  if j.hasKey("collisions"): result.collisions = j["collisions"].getInt
  if j.hasKey("efinal"): result.efinal = j["efinal"].getFloat.eV
  if j.hasKey("penning"): result.penning = j["penning"].getBool
  if j.hasKey("thermal"): result.thermal = j["thermal"].getBool
  if j.hasKey("seed"): result.seed = j["seed"].getInt
  if j.hasKey("report") and j["report"].getBool: result.output = stdout

var answer: string  # lives until the next call, for the caller to read

proc nimboltzRun(json: cstring): cstring {.exportc, cdecl.} =
  ## runs a setup (JSON) and returns its results (JSON)
  try:
    answer = $simulate(setup(parseJson($json))).toJson
  except CatchableError as e:
    answer = $(%*{"error": e.msg})
  stdout.flushFile
  result = cstring(answer)

proc nimboltzGases(): cstring {.exportc, cdecl.} =
  ## the gases, [{"id", "name"}, ...], by Magboltz's gas number
  var a = newJArray()
  for g in Gas:
    a.add %*{"id": ord(g), "name": $g}
  answer = $a
  result = cstring(answer)

proc nimboltzCollisionScale(): cint {.exportc, cdecl.} =
  ## collisions per unit of Magboltz's NMAX (4·10⁷ unless built otherwise)
  cint(nscale)

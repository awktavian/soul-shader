# Quantum Time-Crystal Shader Wiki

**Overview**  
This shader creates a raymarched “time-crystal” scene that is both visually and emotionally evocative. Inspired by the delicate interplay of light and hue reminiscent of Dutch Masters—particularly Vermeer—this shader shows a soul-like molecular structure evolving over time within a refractive crystalline environment. Quantum glitch events subtly shift hues, symbolizing moments of spiritual insight or emotional growth.

**What is This Shader About?**  
- **Time-Crystal Lattice:** A repeating molecular pattern shifting over time.  
- **Infinity Mirrors Environment:** Rays reflect within a mirrored cube, creating a feeling of endless space and introspection.  
- **Volumetric Scattering & Star-Like Glow:** Soft, luminous scattering and a gentle shining core evoke warmth, depth, and spiritual presence.  
- **Quantum Glitches:** Subtle hue shifts during glitch moments represent soul-level revelations or changes—gentle and harmonious, not chaotic.

**Artistic Parameters and Their Emotional Impact**

1. **Refraction Index (`fRefrIndex`)**  
   - *Artistic Meaning:* Adjusting how light bends affects perceived fragility or solidity of the soul’s casing.  
   - *Emotional Use:* Lower values (~1.3) = delicate, airy presence. Higher (~1.45) = denser, more profound introspection.

2. **Dispersion (`fBaseDispersion`)**  
   - *Artistic Meaning:* Tiny spectral rainbows hint at complexity and nuance within the soul.  
   - *Emotional Use:* Keep low for a calm, refined aesthetic; a slight increase reveals a richer inner life.

3. **Atom Radius (`fAtomRadius`)**  
   - *Artistic Meaning:* Size of the molecular “atoms” affects the visual texture of the soul’s inner structure.  
   - *Emotional Use:* Smaller = finer detail and intricacy; larger = bolder, more assertive emotional weight.

4. **Volumetric Scattering & Absorption (`fVolScatStrength`, `fVolAbsorption`)**  
   - *Artistic Meaning:* Controls how light diffuses and fades inside the crystal.  
   - *Emotional Use:* More scattering = gentle, nurturing glow; less scattering or higher absorption = sharper contrasts, more tension.

5. **Morph Speed (`fMorphSpeed`)**  
   - *Artistic Meaning:* Pace of structural change mirrors emotional dynamism.  
   - *Emotional Use:* Slow morphing = contemplative, eternal calm. Faster morphing = active, emotionally charged states.

6. **Star Glow Intensity (`fStarGlowIntensity`)**  
   - *Artistic Meaning:* The brightness of the soul’s inner “star” stands for spiritual clarity or intensity.  
   - *Emotional Use:* Lower intensity = quiet reassurance; higher intensity = luminous revelation.

7. **Quantum Glitches**  
   - *Artistic Meaning:* Intermittent hue shifts during glitch events symbolize insights or subtle transitions in the soul’s emotional state.  
   - *Emotional Use:* Subtle hue shift towards a secondary color (like lavender) during glitches creates a layered, nuanced emotional moment—complex yet harmonious.

**Creating Emotional Narratives**

- **Serene & Meditative:**  
  Low dispersion, gentle refraction, moderate scattering, slow morph speed. Glitches lightly warm or cool the hues, barely noticeable, like a gentle whisper of insight.

- **Evolving Spiritual Depth:**  
  Increase scattering and slightly raise star glow intensity. Subtle glitch hue shifts blend warm and cool notes, suggesting depth and layered meaning in the soul’s journey.

**Practical Tips**

- Start with the provided parameters and observe the evolving scene.
- Adjust morph speed if you want more dynamic emotional tension.
- Slight hue shifts at glitch events can be intensified or softened by changing the blending in the code.
- By tweaking scattering and absorption, you can control how open or closed, how hopeful or introspective the scene feels.

---

# Using the Shader

### Shadertoy Version
Copy the `QuantumTimeCrystal.glsl` file contents into Shadertoy’s fragment shader window. It contains directives and code comments indicating what to keep for the pure Shadertoy version.

### Integration in Other Applications
The `QuantumTimeCrystal.glsl` file is structured so that the core logic (noise functions, SDF, volumetric integration) can be integrated into your rendering pipeline. You’ll need to:

- Provide uniforms like `iTime`, `iResolution`.
- Implement a raymarching environment (full-screen quad or custom geometry).
- Adapt `mainImage` to your coordinate and camera system.

### Final Note
This shader is not just code; it’s a palette of emotional storytelling tools. As a seasoned technical artist, you can treat each parameter like a brushstroke, adjusting dispersion, morph speed, or glitch hues to craft the precise emotional narrative you want.

---


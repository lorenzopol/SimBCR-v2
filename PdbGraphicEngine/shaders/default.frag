#version 330

out vec4 finalColor;

in vec2 fragTexCoord;
in vec3 fragNormal;
in vec3 fragPosition;


uniform sampler2D texture0;
uniform vec3 viewPos;

uniform vec3 ambient_intensity;
uniform vec3 diffuse_intensity;
uniform vec3 specular_intensity;
uniform vec3 light_pos;


vec3 getLight(vec3 color) {
    vec3 Normal = normalize(fragNormal);

    // diffuse light (Lambert law)
    vec3 lightDir = normalize(light_pos - fragPosition);
    float diff  = max(0, dot(lightDir, Normal));
    vec3 diffuse = diff * diffuse_intensity;

    // specular
    vec3 viewDir = normalize(viewPos - fragPosition);
    vec3 reflectDir = reflect(-lightDir, Normal);
    // values bigger than 32 will make the surface shinier, less then 32 will make it roughier
    float spec = pow(max(dot(viewDir, reflectDir), 0), 32);
    vec3 specular = spec * specular_intensity;

    return color * (ambient_intensity + diffuse + specular);
}
void main() {
    float gamma = 2.2;

    vec3 start_color = vec3(1.0, 0.0, 0.0);
    //vec3 start_color = texture(texture0, fragTexCoord).rgb;
    // 1. undo gamma correction (done in the tex)
    vec3 color = pow(start_color, vec3(gamma));

    // 2. calculate light
    color = getLight(color);

    // 3. gamma correct the blended light + tex
    color = pow(color, 1/vec3(gamma));
    finalColor = vec4(color, 1.0);
}
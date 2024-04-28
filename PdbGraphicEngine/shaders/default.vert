#version 330

in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec3 vertexPosition;
in mat4 instanceTransform;

uniform mat4 mvp;
uniform mat4 matNormal;

out vec2 fragTexCoord;
out vec3 fragNormal;
out vec3 fragPosition;


void main() {
    fragTexCoord = vertexTexCoord;
    mat4 mvpi = mvp*instanceTransform;

    // position of the fragment in world space
    fragPosition = vec3(mvpi * vec4(vertexPosition, 1.0));

    // transpose(invert(m_model)) is done in order to prevent normal to be affectaed by non uniform scale
    fragNormal = mat3(transpose(inverse(mvpi))) * normalize(vertexNormal);
    gl_Position = mvpi * vec4(vertexPosition, 1.0);
}